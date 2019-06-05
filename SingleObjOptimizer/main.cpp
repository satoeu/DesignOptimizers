#define _CRT_SECURE_NO_WARNINGS

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <omp.h>
#include <fstream>
#include <string>
using namespace std;

#ifdef MPI_RUN_OPT_MODE_SETTINGS
#include <mpi.h>
#endif

#include <BaseDefines.hpp>
#include <BaseFuncs/Mt.hpp>

#include "01_OptSetting.hpp"
#include "RPopulationS.hpp"


#include <TestObjFuncs/TestObjFuncs.hpp>

#include "SingleObjOptimizer.hpp"


#ifdef IS_WINDOWS_SISTEM
#ifdef _DEBUG
#pragma comment(lib, "libBaseFuncs_Deb.lib")
#pragma comment(lib, "TestObjFuncs_Deb.lib")
#else
#pragma comment(lib, "libBaseFuncs.lib")
#pragma comment(lib, "TestObjFuncs.lib")
#endif
#endif




/*//=======================================================
  // ●　メイン関数
  //=======================================================*/
int main(int argc, char *argv[]){

#ifdef MPI_RUN_OPT_MODE_SETTINGS
	/* MPI初期化 */
	MPI_Init(&argc, &argv);
	int mpi_rank, mpi_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif



	int seed = (unsigned)time(NULL);

#ifdef MPI_RUN_OPT_MODE_SETTINGS
	if(mpi_rank == 0){
		cout << "Seed: "<< seed << endl;
	}
#endif

	CommonLibs::Mt::init_rand(seed);
	srand(seed);

	/* 設定初期化 */
	OptSettings::opt_set_read();
	

	/* 実行 */
	SingleObjOptimizer optimizer;

	
#ifdef MPI_RUN_OPT_MODE_SETTINGS
	optimizer.run_optimization_MPI();
#else
	optimizer.run_optimization();
#endif


#ifdef MPI_RUN_OPT_MODE_SETTINGS
	/* MPI終了 */
	MPI_Finalize();
#else
	return 1;
#endif
}



