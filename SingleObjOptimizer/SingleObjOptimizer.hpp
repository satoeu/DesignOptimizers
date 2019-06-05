#ifndef DEF_OPT_DIRVER
#define DEF_OPT_DIRVER

#include <omp.h>
#include <BaseDefines.hpp>
#include "01_OptSetting.hpp"
#include "RPopulationS.hpp"
#include "FunctionalOracle.hpp"

#ifdef MPI_RUN_OPT_MODE_SETTINGS
#include <mpi.h>
#endif

#include "Generators/AREX.hpp"
#include "Generators/REXstar.hpp"
#include "Optimizers/JGG.hpp"

/*
//=============================================================
// ■ SingleObjOptimizer
//=============================================================
//	単目的最適化実行主体クラス
//============================================================= */
class SingleObjOptimizer{
private:
    RPopulationS population;								/* 個体集団 */
    FunctionalOracle* oracles;                          /* 評価実施用関数オラクル */
    int total_eval_count;                               /* 総評価回数セーブ */
    /**/
    void eval_new_genes(vector<RGeneS>& rgenes);         /* 個体群の評価を実施評価 */
	//void init_oracle_ref(const int para_num, double*** gene_values, int** conv_flag, double*** oracl_results);		/* Oracle実行に必要な情報の初期化 */
	//void del_oracle_ref(const int para_num, double*** gene_values, int** conv_flag, double*** oracl_results);		/* Oracle実行に必要な情報の削除 */
#ifdef MPI_RUN_OPT_MODE_SETTINGS
    void eval_new_genes_MPI(const int eval_size, vector<RGeneS>& rgenes);								/* 個体群の評価を実施評価(MPIのとき) */
	void calc_init_population_MPIver();																   /* 初期個体群を評価(MPIver)  */
#endif
/**/
public:
	SingleObjOptimizer();                               /* コンストラクタ */
    ~SingleObjOptimizer();                              /* デストラクタ */
    void run_optimization();							/* 最適化実行 */
#ifdef MPI_RUN_OPT_MODE_SETTINGS
    void run_optimization_MPI();                        /* 最適化実行(MPIのとき) */
#endif
};


#endif
