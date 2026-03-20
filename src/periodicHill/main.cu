// ================================================================
// Periodic Hill Flow — Main Simulation Program (FluidX3D Edition)
// 對齊 GILBM D3Q19 reference code 的 main.cu 結構
//
// 此檔案被 setup.cpp #include，並定義 main_setup() 函數
// 所有 FluidX3D API (LBM, Memory_Container, utilities) 已透過 setup.hpp 導入
// ================================================================

#include <iomanip>
#include <cstdio>
#include <fstream>
#include <cstring>
#include <cmath>

// ================================================================
// Hill Geometry (Mellen et al. piecewise cubic polynomial)
// ================================================================
static float hill_poly(const float s) {
	const float s2=s*s, s3=s2*s;
	if(s<=9.0f) {
		return fmin(28.0f, 28.0f+0.006775070969851f*s2-0.0021245277758f*s3);
	} else if(s<=14.0f) {
		return 25.07355893131f+0.9754803562315f*s-0.1016116352781f*s2+0.001889794677828f*s3;
	} else if(s<=20.0f) {
		return 25.79601052357f+0.8206693007457f*s-0.09055370274339f*s2+0.001626510569859f*s3;
	} else if(s<=30.0f) {
		return 40.46435022819f-1.379581654948f*s+0.019458845041284f*s2-0.000207031893219f*s3;
	} else if(s<=40.0f) {
		return 17.92461334664f+0.8743920332081f*s-0.05567361123058f*s2+0.0006277731764683f*s3;
	} else if(s<=54.0f) {
		return fmax(0.0f, 56.39011190988f-2.010520359035f*s+0.01644919857549f*s2+0.00002674976141766f*s3);
	}
	return 0.0f;
}

static float hill_height(const float y, const float Ny, const float h) {
	const float LY = (float)PH_LY;
	const float Y = y*LY/Ny;
	const float s_left  = Y*28.0f;
	const float s_right = (LY-Y)*28.0f;
	float poly = 0.0f;
	if(s_left<=54.0f) {
		poly = hill_poly(s_left);
	} else if(s_right<=54.0f) {
		poly = hill_poly(s_right);
	}
	return poly/28.0f*h;
}

// ================================================================
// Binary VTK Writer (Big-Endian, aligned with FluidX3D native format)
// ================================================================
// Writes ./result/velocity_merged_XXXXXX.vtk
// Before FTT_STATS_START: instantaneous fields only
// After  FTT_STATS_START: instantaneous + mean + Reynolds stresses
static void write_vtk_binary(
	LBM& lbm, const ulong step, const float Uref,
	const ulong accu_count,
	const double* sum_ux, const double* sum_uy, const double* sum_uz,
	const double* sum_uy2, const double* sum_uz2,
	const double* sum_ux2, const double* sum_uyuz
) {
	const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz();
	const ulong N = lbm.get_N();

	char fname[256];
	snprintf(fname, sizeof(fname), "./result/velocity_merged_%06llu.vtk", (unsigned long long)step);

	// ASCII header
	std::string header;
	{
		char buf[512];
		snprintf(buf, sizeof(buf),
			"# vtk DataFile Version 3.0\n"
			"PeriodicHill step=%llu Force=%.8e accu_count=%llu\n"
			"BINARY\n"
			"DATASET STRUCTURED_POINTS\n"
			"DIMENSIONS %u %u %u\n"
			"ORIGIN 0 0 0\n"
			"SPACING 1 1 1\n"
			"POINT_DATA %llu\n",
			(unsigned long long)step, (double)lbm.get_fy(), (unsigned long long)accu_count,
			Nx, Ny, Nz, (unsigned long long)N);
		header = buf;
	}

	create_folder("./result/dummy.txt");
	std::ofstream file(fname, std::ios::out|std::ios::binary);
	if(!file.is_open()) { print_info("ERROR: Cannot open VTK file: "+string(fname)); return; }
	file.write(header.c_str(), header.length());

	// helper: write a SCALARS block (binary big-endian float)
	auto write_scalar = [&](const char* name, const float* data_f, const ulong count) {
		char hdr[128];
		snprintf(hdr, sizeof(hdr), "SCALARS %s float 1\nLOOKUP_TABLE default\n", name);
		file.write(hdr, strlen(hdr));
		for(ulong i=0; i<count; i++) {
			float val = reverse_bytes(data_f[i]);
			file.write((char*)&val, sizeof(float));
		}
	};

	// --- VECTORS velocity (benchmark coords: streamwise=y, spanwise=x, wall-normal=z) ---
	{
		char hdr[] = "VECTORS velocity float\n";
		file.write(hdr, strlen(hdr));
		const float inv_Uref = 1.0f/Uref;
		for(ulong n=0; n<N; n++) {
			float v[3] = {
				reverse_bytes(lbm.u.y[n]*inv_Uref),  // streamwise (benchmark U)
				reverse_bytes(lbm.u.x[n]*inv_Uref),  // spanwise   (benchmark V)
				reverse_bytes(lbm.u.z[n]*inv_Uref)    // wall-normal(benchmark W)
			};
			file.write((char*)v, 3*sizeof(float));
		}
	}

	// --- Instantaneous velocity scalars ---
	{
		const float inv_Uref = 1.0f/Uref;
		float* tmp = new float[N];

		for(ulong n=0; n<N; n++) tmp[n] = lbm.u.y[n]*inv_Uref; // u_inst = streamwise
		write_scalar("u_inst", tmp, N);

		for(ulong n=0; n<N; n++) tmp[n] = lbm.u.x[n]*inv_Uref; // v_inst = spanwise
		write_scalar("v_inst", tmp, N);

		for(ulong n=0; n<N; n++) tmp[n] = lbm.u.z[n]*inv_Uref; // w_inst = wall-normal
		write_scalar("w_inst", tmp, N);

		delete[] tmp;
	}

	// --- Time-averaged fields (only after FTT_STATS_START, when accu_count > 0) ---
	if(accu_count>0 && sum_uy!=nullptr) {
		const double inv_count = 1.0/(double)accu_count;
		const double inv_Uref2 = 1.0/((double)Uref*(double)Uref);
		float* tmp = new float[N];

		// U_mean (streamwise mean / Uref)
		for(ulong n=0; n<N; n++) tmp[n] = (float)(sum_uy[n]*inv_count/(double)Uref);
		write_scalar("U_mean", tmp, N);

		// W_mean (wall-normal mean / Uref)
		for(ulong n=0; n<N; n++) tmp[n] = (float)(sum_uz[n]*inv_count/(double)Uref);
		write_scalar("W_mean", tmp, N);

#if PH_VTK_OUTPUT_LEVEL >= 1
		// V_mean (spanwise mean / Uref)
		for(ulong n=0; n<N; n++) tmp[n] = (float)(sum_ux[n]*inv_count/(double)Uref);
		write_scalar("V_mean", tmp, N);
#endif

		// uu_RS = (<uy²> - <uy>²) / Uref²
		for(ulong n=0; n<N; n++) {
			double mean_uy = sum_uy[n]*inv_count;
			tmp[n] = (float)((sum_uy2[n]*inv_count - mean_uy*mean_uy)*inv_Uref2);
		}
		write_scalar("uu_RS", tmp, N);

		// uw_RS = (<uy·uz> - <uy><uz>) / Uref²
		for(ulong n=0; n<N; n++) {
			double mean_uy = sum_uy[n]*inv_count, mean_uz = sum_uz[n]*inv_count;
			tmp[n] = (float)((sum_uyuz[n]*inv_count - mean_uy*mean_uz)*inv_Uref2);
		}
		write_scalar("uw_RS", tmp, N);

		// ww_RS = (<uz²> - <uz>²) / Uref²
		for(ulong n=0; n<N; n++) {
			double mean_uz = sum_uz[n]*inv_count;
			tmp[n] = (float)((sum_uz2[n]*inv_count - mean_uz*mean_uz)*inv_Uref2);
		}
		write_scalar("ww_RS", tmp, N);

#if PH_VTK_OUTPUT_LEVEL >= 1
		// vv_RS (spanwise) = (<ux²> - <ux>²) / Uref²
		for(ulong n=0; n<N; n++) {
			double mean_ux = sum_ux[n]*inv_count;
			tmp[n] = (float)((sum_ux2[n]*inv_count - mean_ux*mean_ux)*inv_Uref2);
		}
		write_scalar("vv_RS", tmp, N);
#endif

		// k_TKE = 0.5*(uu + vv + ww) / Uref²
		for(ulong n=0; n<N; n++) {
			double mux = sum_ux[n]*inv_count, muy = sum_uy[n]*inv_count, muz = sum_uz[n]*inv_count;
			double uu=(sum_uy2[n]*inv_count-muy*muy)*inv_Uref2;
			double vv=(sum_ux2[n]*inv_count-mux*mux)*inv_Uref2;
			double ww=(sum_uz2[n]*inv_count-muz*muz)*inv_Uref2;
			tmp[n] = (float)(0.5*(uu+vv+ww));
		}
		write_scalar("k_TKE", tmp, N);

		delete[] tmp;
	}

	file.close();
	print_info("VTK binary: "+string(fname)+(accu_count>0?" (accu="+to_string(accu_count)+")":""));
}

// ================================================================
// Binary Checkpoint I/O
// ================================================================
// 目錄結構: checkpoint/step_XXXXXX/
//   metadata.dat  — step, FTT, accu_count, Force, tau
//   rho.bin       — density           [N × float]
//   ux.bin        — velocity x        [N × float]
//   uy.bin        — velocity y        [N × float]
//   uz.bin        — velocity z        [N × float]
//   flags.bin     — cell flags        [N × uchar]
//   (Stage 1 only, accu_count > 0):
//   sum_ux.bin ... sum_uyuz.bin — 統計量累積和 [N × double]
//   accu_count.dat — 累積次數
//
// 注意: FluidX3D 的分佈函數 (DDFs) 為 device-only memory，
//       無法從 setup.cpp 存取。因此 checkpoint 儲存巨觀量 (rho, u)，
//       重啟時使用 feq(rho, u) 近似，等同 reference code 的 INIT=2 模式。
// ================================================================

static void write_binary_array_float(const string& filepath, const float* data, const ulong N) {
	std::ofstream f(filepath, std::ios::binary);
	f.write((const char*)data, N*sizeof(float));
	f.close();
}

static void write_binary_array_double(const string& filepath, const double* data, const ulong N) {
	std::ofstream f(filepath, std::ios::binary);
	f.write((const char*)data, N*sizeof(double));
	f.close();
}

static void write_binary_array_uchar(const string& filepath, const uchar* data, const ulong N) {
	std::ofstream f(filepath, std::ios::binary);
	f.write((const char*)data, N*sizeof(uchar));
	f.close();
}

static bool read_binary_array_float(const string& filepath, float* data, const ulong N) {
	std::ifstream f(filepath, std::ios::binary);
	if(!f.is_open()) return false;
	f.read((char*)data, N*sizeof(float));
	f.close();
	return true;
}

static bool read_binary_array_double(const string& filepath, double* data, const ulong N) {
	std::ifstream f(filepath, std::ios::binary);
	if(!f.is_open()) return false;
	f.read((char*)data, N*sizeof(double));
	f.close();
	return true;
}

static void save_checkpoint(
	LBM& lbm, const ulong step, const double FTT,
	const ulong accu_count,
	const double* sum_ux, const double* sum_uy, const double* sum_uz,
	const double* sum_uy2, const double* sum_uz2,
	const double* sum_ux2, const double* sum_uyuz
) {
	const ulong N = lbm.get_N();
	char dir[256];
	snprintf(dir, sizeof(dir), "checkpoint/step_%06llu", (unsigned long long)step);

	// Create directories
	create_folder("checkpoint/dummy.txt");
	create_folder(string(dir)+"/dummy.txt");

	// Write macroscopic fields (equivalent to distribution function checkpoint for feq restart)
	// rho
	{
		float* buf = new float[N];
		for(ulong n=0; n<N; n++) buf[n] = lbm.rho[n];
		write_binary_array_float(string(dir)+"/rho.bin", buf, N);
		delete[] buf;
	}
	// ux, uy, uz
	{
		float* buf = new float[N];
		for(ulong n=0; n<N; n++) buf[n] = lbm.u.x[n];
		write_binary_array_float(string(dir)+"/ux.bin", buf, N);
		for(ulong n=0; n<N; n++) buf[n] = lbm.u.y[n];
		write_binary_array_float(string(dir)+"/uy.bin", buf, N);
		for(ulong n=0; n<N; n++) buf[n] = lbm.u.z[n];
		write_binary_array_float(string(dir)+"/uz.bin", buf, N);
		delete[] buf;
	}
	// flags
	{
		uchar* buf = new uchar[N];
		for(ulong n=0; n<N; n++) buf[n] = lbm.flags[n];
		write_binary_array_uchar(string(dir)+"/flags.bin", buf, N);
		delete[] buf;
	}

	// Statistics arrays (only if accumulation has started)
	if(accu_count > 0 && sum_uy != nullptr) {
		write_binary_array_double(string(dir)+"/sum_ux.bin",   sum_ux,   N);
		write_binary_array_double(string(dir)+"/sum_uy.bin",   sum_uy,   N);
		write_binary_array_double(string(dir)+"/sum_uz.bin",   sum_uz,   N);
		write_binary_array_double(string(dir)+"/sum_uy2.bin",  sum_uy2,  N);
		write_binary_array_double(string(dir)+"/sum_uz2.bin",  sum_uz2,  N);
		write_binary_array_double(string(dir)+"/sum_ux2.bin",  sum_ux2,  N);
		write_binary_array_double(string(dir)+"/sum_uyuz.bin", sum_uyuz, N);
	}

	// metadata.dat
	{
		char meta_path[256];
		snprintf(meta_path, sizeof(meta_path), "%s/metadata.dat", dir);
		std::ofstream meta(meta_path);
		meta << "step=" << step << "\n";
		meta << std::fixed << std::setprecision(15);
		meta << "FTT=" << FTT << "\n";
		meta << "accu_count=" << accu_count << "\n";
		meta << "Force=" << (double)lbm.get_fy() << "\n";
		meta << "tau=" << (double)lbm.get_tau() << "\n";
		meta << "Nx=" << lbm.get_Nx() << "\n";
		meta << "Ny=" << lbm.get_Ny() << "\n";
		meta << "Nz=" << lbm.get_Nz() << "\n";
		meta.close();
	}

	print_info("[CHECKPOINT] Saved: "+string(dir)+"/ (FTT="+to_string((float)FTT, 2u)+", accu="+to_string(accu_count)+")");
}

static bool load_checkpoint(
	LBM& lbm, const char* checkpoint_dir,
	ulong& out_accu_count,
	double* sum_ux, double* sum_uy, double* sum_uz,
	double* sum_uy2, double* sum_uz2,
	double* sum_ux2, double* sum_uyuz
) {
	const ulong N = lbm.get_N();
	string dir(checkpoint_dir);

	// Read metadata
	{
		std::ifstream meta(dir+"/metadata.dat");
		if(!meta.is_open()) {
			print_info("ERROR: Cannot open "+dir+"/metadata.dat");
			return false;
		}
		string line;
		while(std::getline(meta, line)) {
			if(line.find("accu_count=")!=string::npos) {
				out_accu_count = (ulong)std::stoull(line.substr(line.find('=')+1));
			}
		}
		meta.close();
	}

	// Read macroscopic fields → write back into lbm host arrays
	{
		float* buf = new float[N];
		if(read_binary_array_float(dir+"/rho.bin", buf, N)) {
			for(ulong n=0; n<N; n++) lbm.rho[n] = buf[n];
		}
		if(read_binary_array_float(dir+"/ux.bin", buf, N)) {
			for(ulong n=0; n<N; n++) lbm.u.x[n] = buf[n];
		}
		if(read_binary_array_float(dir+"/uy.bin", buf, N)) {
			for(ulong n=0; n<N; n++) lbm.u.y[n] = buf[n];
		}
		if(read_binary_array_float(dir+"/uz.bin", buf, N)) {
			for(ulong n=0; n<N; n++) lbm.u.z[n] = buf[n];
		}
		delete[] buf;
	}

	// Read statistics (if file exists)
	if(out_accu_count > 0) {
		read_binary_array_double(dir+"/sum_ux.bin",   sum_ux,   N);
		read_binary_array_double(dir+"/sum_uy.bin",   sum_uy,   N);
		read_binary_array_double(dir+"/sum_uz.bin",   sum_uz,   N);
		read_binary_array_double(dir+"/sum_uy2.bin",  sum_uy2,  N);
		read_binary_array_double(dir+"/sum_uz2.bin",  sum_uz2,  N);
		read_binary_array_double(dir+"/sum_ux2.bin",  sum_ux2,  N);
		read_binary_array_double(dir+"/sum_uyuz.bin", sum_uyuz, N);
	}

	print_info("[CHECKPOINT] Loaded: "+dir+"/ (accu="+to_string(out_accu_count)+")");
	return true;
}

// ================================================================
// main_setup() — Periodic Hill Flow 專案主程式
// ================================================================
void main_setup() {
	// ##################### 1. 建立 LBM 物件 #####################
	LBM lbm((uint)PH_NX, (uint)PH_NY, (uint)PH_NZ, PH_niu, 0.0f, PH_FY, 0.0f);

	const uint Nx_ = lbm.get_Nx();
	const uint Ny_ = lbm.get_Ny();
	const uint Nz_ = lbm.get_Nz();
	const ulong N  = lbm.get_N();

	// ##################### 2. 定義幾何 #####################
	parallel_for(N, [&](ulong n) {
		uint x=0u, y=0u, z=0u;
		lbm.coordinates(n, x, y, z);
		const float hz = hill_height((float)y, (float)Ny_, (float)PH_H);
		if((float)z <= hz || z >= Nz_-1u) {
			lbm.flags[n] = TYPE_S;
		} else {
			lbm.u.y[n] = PH_Uref; // 初始化流向速度，加速收斂
		}
	});

	// ##################### 3. 顯示模擬參數 #####################
	const float tau = 3.0f*PH_niu + 0.5f;
	const float Ma  = PH_Uref / 0.57735f;
	const double flow_through_time = (double)Ny_ / (double)PH_Uref;
	const ulong total_steps = (ulong)(PH_FTT_STOP * flow_through_time);

	print_info("+================================================================+");
	print_info("| Periodic Hill Flow — FluidX3D                                  |");
	print_info("+================================================================+");
	print_info("| Grid:  Nx="+to_string(Nx_)+", Ny="+to_string(Ny_)+", Nz="+to_string(Nz_)+" (N="+to_string(N)+")");
	print_info("| Hill:  h="+to_string(PH_H)+" cells, Re="+to_string(PH_Re));
	print_info("| Uref="+to_string(PH_Uref, 6u)+", nu="+to_string(PH_niu, 6u)+", fy="+to_string(PH_FY, 8u));
	print_info("| tau="+to_string(tau, 4u)+", Ma="+to_string(Ma, 4u));
	print_info("| FTT="+to_string((float)flow_through_time, 1u)+" steps, total="+to_string(total_steps)+" steps");
	print_info("| NDTMIT="+to_string(PH_NDTMIT)+", NDTVTK="+to_string(PH_NDTVTK)+", NDTCKPT="+to_string(PH_NDTCKPT));
	print_info("| FTT_STATS_START="+to_string((float)PH_FTT_STATS_START, 1u)+", FTT_STOP="+to_string((float)PH_FTT_STOP, 1u));
	print_info("+================================================================+");

#ifdef GRAPHICS
	lbm.graphics.visualization_modes = VIS_FLAG_LATTICE|VIS_FIELD;
	lbm.graphics.slice_mode = 2;
	lbm.graphics.slice_y = (int)(Ny_/2u);
#endif

	// ##################### 4. 分配統計量陣列 (host, double) #####################
	double* sum_ux   = new double[N]();
	double* sum_uy   = new double[N]();
	double* sum_uz   = new double[N]();
	double* sum_uy2  = new double[N]();
	double* sum_uz2  = new double[N]();
	double* sum_ux2  = new double[N]();
	double* sum_uyuz = new double[N]();
	ulong accu_count = 0u;

	// ##################### 5. Restart (if INIT=1) #####################
#if PH_INIT == 1
	{
		print_info("Loading checkpoint: " PH_RESTART_DIR);
		bool ok = load_checkpoint(lbm, PH_RESTART_DIR, accu_count,
			sum_ux, sum_uy, sum_uz, sum_uy2, sum_uz2, sum_ux2, sum_uyuz);
		if(!ok) {
			print_info("WARNING: Checkpoint load failed, starting from scratch.");
			accu_count = 0;
		} else {
			// Write loaded velocity back to device (feq restart)
			lbm.rho.write_to_device();
			lbm.u.write_to_device();
		}
	}
#endif

	// ##################### 6. 初始化 monitor 檔案 #####################
	create_folder("./result/dummy.txt");
	create_folder("./checkpoint/dummy.txt");
	{
		std::ofstream f("checkrho.dat", std::ios::trunc);
		f << "# step\tFTT\trho_ref\trho_mean\n";
		f.close();
	}
	{
		std::ofstream f("Ustar_Force_record.dat", std::ios::trunc);
		f << "# FTT\tUb/Uref\tForce\tMa_max\taccu_count\n";
		f.close();
	}
	// Log file
	{
		std::ofstream f("checkpoint/log.dat", std::ios::trunc);
		f << "# Periodic Hill Flow — FluidX3D Run Log\n";
		f << "# step\tFTT\tUb/Uref\tMa_max\trho_mean\taccu_count\n";
		f.close();
	}

	// ##################### 7. 主時間迴圈 #####################
	while(true) {
		lbm.run((ulong)PH_NDTMIT, total_steps);
		const ulong step = lbm.get_t();
		const double FTT_now = (double)step / flow_through_time;

		// --- GPU → Host 資料傳輸 ---
		lbm.u.read_from_device();
		lbm.rho.read_from_device();
		lbm.flags.read_from_device();

		// --- 計算全場 bulk velocity (流向平均) ---
		double Ub_sum = 0.0;
		ulong Ub_count = 0u;
		for(ulong n=0; n<N; n++) {
			if(lbm.flags[n] != TYPE_S) {
				Ub_sum += (double)lbm.u.y[n];
				Ub_count++;
			}
		}
		const double Ub_inst = (Ub_count>0) ? Ub_sum/(double)Ub_count : 0.0;

		// --- 計算最大 Mach 數 ---
		float max_u_sq = 0.0f;
		for(ulong n=0; n<N; n++) {
			const float usq = sq(lbm.u.x[n]) + sq(lbm.u.y[n]) + sq(lbm.u.z[n]);
			if(usq > max_u_sq) max_u_sq = usq;
		}
		const double Ma_max = (double)sqrt(max_u_sq) / 0.57735;

		// --- 計算平均密度 ---
		double rho_sum = 0.0;
		ulong rho_count = 0u;
		for(ulong n=0; n<N; n++) {
			if(lbm.flags[n] != TYPE_S) {
				rho_sum += (double)lbm.rho[n];
				rho_count++;
			}
		}
		const double rho_mean = (rho_count>0) ? rho_sum/(double)rho_count : 1.0;

		// --- checkrho.dat (每 NDTMIT 步) ---
		{
			std::ofstream f("checkrho.dat", std::ios::app);
			f << step << "\t" << std::fixed << std::setprecision(4) << FTT_now
			  << "\t" << std::setprecision(6) << 1.0 << "\t" << rho_mean << "\n";
			f.close();
		}

		// --- Ustar_Force_record.dat (每 NDTMIT 步) ---
		{
			const double F_star = (double)lbm.get_fy() * (double)Ny_ / ((double)PH_Uref*(double)PH_Uref);
			std::ofstream f("Ustar_Force_record.dat", std::ios::app);
			f << std::fixed << std::setprecision(6) << FTT_now << "\t"
			  << std::setprecision(10) << Ub_inst/(double)PH_Uref << "\t"
			  << std::setprecision(10) << F_star << "\t"
			  << std::setprecision(6)  << Ma_max << "\t"
			  << accu_count << "\n";
			f.close();
		}

		// --- Log 紀錄 (每 NDTMIT 步) ---
		{
			std::ofstream f("checkpoint/log.dat", std::ios::app);
			f << step << "\t"
			  << std::fixed << std::setprecision(4) << FTT_now << "\t"
			  << std::setprecision(6) << Ub_inst/(double)PH_Uref << "\t"
			  << std::setprecision(4) << Ma_max << "\t"
			  << std::setprecision(6) << rho_mean << "\t"
			  << accu_count << "\n";
			f.close();
		}

		// --- 累積統計量 (FTT ≥ FTT_STATS_START) ---
		if(FTT_now >= PH_FTT_STATS_START) {
			for(ulong n=0; n<N; n++) {
				const double ux=(double)lbm.u.x[n], uy=(double)lbm.u.y[n], uz=(double)lbm.u.z[n];
				sum_ux[n]   += ux;
				sum_uy[n]   += uy;
				sum_uz[n]   += uz;
				sum_uy2[n]  += uy*uy;
				sum_uz2[n]  += uz*uz;
				sum_ux2[n]  += ux*ux;
				sum_uyuz[n] += uy*uz;
			}
			accu_count++;
		}

		// --- Console 進度顯示 ---
		if(step % (ulong)PH_PRINT_INTERVAL == 0u) {
			print_info("  step="+to_string(step)
				+", FTT="+to_string((float)FTT_now, 2u)
				+", Ub/Uref="+to_string((float)(Ub_inst/(double)PH_Uref), 6u)
				+", Ma_max="+to_string((float)Ma_max, 4u)
				+", rho="+to_string((float)rho_mean, 6u)
				+(accu_count>0 ? ", accu="+to_string(accu_count) : ""));
		}

		// --- VTK 輸出 (每 NDTVTK 步, binary) ---
		if(step % (ulong)PH_NDTVTK == 0u) {
			write_vtk_binary(lbm, step, PH_Uref, accu_count,
				sum_ux, sum_uy, sum_uz, sum_uy2, sum_uz2, sum_ux2, sum_uyuz);
		}

		// --- Binary checkpoint (每 NDTCKPT 步) ---
		// FTT < FTT_STATS_START: 只寫 rho + u (分佈函數近似)
		// FTT ≥ FTT_STATS_START: 同時寫入統計量累積和
		if(step % (ulong)PH_NDTCKPT == 0u) {
			save_checkpoint(lbm, step, FTT_now, accu_count,
				sum_ux, sum_uy, sum_uz, sum_uy2, sum_uz2, sum_ux2, sum_uyuz);
		}

		// --- FTT 停止判定 ---
		if(FTT_now >= PH_FTT_STOP) {
			print_info("[FTT-STOP] FTT="+to_string((float)FTT_now, 2u)+" >= FTT_STOP="+to_string((float)PH_FTT_STOP, 1u)
				+" at step "+to_string(step)+". Ending simulation.");
			// Final VTK + checkpoint
			write_vtk_binary(lbm, step, PH_Uref, accu_count,
				sum_ux, sum_uy, sum_uz, sum_uy2, sum_uz2, sum_ux2, sum_uyuz);
			save_checkpoint(lbm, step, FTT_now, accu_count,
				sum_ux, sum_uy, sum_uz, sum_uy2, sum_uz2, sum_ux2, sum_uyuz);
			break;
		}
	}

	// ##################### 8. 清理 #####################
	delete[] sum_ux;
	delete[] sum_uy;
	delete[] sum_uz;
	delete[] sum_uy2;
	delete[] sum_uz2;
	delete[] sum_ux2;
	delete[] sum_uyuz;
	print_info("Periodic Hill simulation complete. accu_count="+to_string(accu_count));
}
