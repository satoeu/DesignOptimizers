#include "SingleObjOptimizer.hpp"
#include <fstream>

/*
//=============================================================
// ■ SingleObjOptimizer
//=============================================================
//	最適化実行主体クラス
//============================================================= */

/*//===========================================================
// ● コンストラクタ
//=========================================================== */
SingleObjOptimizer::SingleObjOptimizer(){

	/* OpenMP並列数セット */
#ifndef MPI_RUN_OPT_MODE_SETTINGS
	int max_thread = omp_get_max_threads();
	max_thread = OptSettings::PARA_OMP_CALC_NUM > max_thread ? max_thread : OptSettings::PARA_OMP_CALC_NUM;
	omp_set_num_threads(max_thread);
#else
	omp_set_num_threads(1);
	int mpi_rank, mpi_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif

	/* 集団初期化 */
	population.initialize(OptSettings::TOTAL_SIZE);

	/* オラクル初期化 */
#ifndef MPI_RUN_OPT_MODE_SETTINGS
	const int PARA_OMP_CALC_NUM = OptSettings::PARA_OMP_CALC_NUM;
	oracles = new FunctionalOracle[PARA_OMP_CALC_NUM];
#else
	oracles = new FunctionalOracle;
#endif

	this->total_eval_count = 0;
}

/*//===========================================================
// ● デストラクタ
//=========================================================== */
SingleObjOptimizer::~SingleObjOptimizer(){
#ifndef MPI_RUN_OPT_MODE_SETTINGS
	delete[] oracles;
#else
	delete oracles;
#endif
}

/*//===========================================================
// ● 最適化実行
//=========================================================== */
void SingleObjOptimizer::run_optimization(){
	/*-------------------------------------------------------*/
	/* 定数群 */
	const int TOTAL_SIZE = OptSettings::TOTAL_SIZE;
	const int GENE_SIZE = OptSettings::GENE_SIZE;
	const int RGENE_FEATURES = OptSettings::get_total_oracle_features();	
	/*-------------------------------------------------------*/
	/* 定数セット */
	const int PARA_OMP_CALC_NUM = OptSettings::PARA_OMP_CALC_NUM;
	const int GENERATION_SIZE = OptSettings::GENERATION_SIZE;
	/*-------------------------------------------------------*/

	fstream fpX;
	fpX.open("GA_adapt.csv", std::ios::out);
	fpX << "generation, num_eval, adapt, ";
	for(int k=0 ; k < OptSettings::NUM_FEATURE ; k++){
		fpX << OptSettings::FEATURES_NAME[k] << ", ";
	}
	fpX << "residual" << endl;
	fpX.close();

	/* OpenMP並列数セット */
	int max_thread = omp_get_max_threads();
	max_thread = PARA_OMP_CALC_NUM > max_thread ? max_thread : PARA_OMP_CALC_NUM;
	omp_set_num_threads(max_thread);

	/* 初期個体群評価 */
	this->eval_new_genes(population.genes);

	/*-------------------------------------------------------*/
#if defined(USING_CROSS_REX_STAR) || defined(USING_CROSS_AREX)
	OptSettings::REX_Settings::opt_set_read_rex();
#endif
#ifdef USING_CROSS_AREX
	SimpleAREX rex_cross(OptSettings::REX_Settings::NUM_CHILD, OptSettings::REX_Settings::PARENT_NUM);
#elif defined(USING_CROSS_REX_STAR)
	REXstar rex_cross(OptSettings::REX_Settings::NUM_CHILD, OptSettings::REX_Settings::PARENT_NUM, 0.1);
#endif
	/*-------------------------------------------------------*/
#if defined(USING_CROSS_REX_STAR) || defined(USING_CROSS_AREX)
	/* 世代交代 */
	JGG jgg(TOTAL_SIZE);
#endif
	/*-------------------------------------------------------*/

	/* 世代ループ */
	for(int gg=0 ; gg < GENERATION_SIZE ; gg++){
		cout << "generation " << gg << endl;
		/*-------------------------------------------------------*/
		/* REX */
		rex_cross.create_children(population.genes, population.children);
		/* 子個体評価 */
		this->eval_new_genes(population.children);
		/*-------------------------------------------------------*/
		/* JGGで個体群更新 */
		jgg.selection(rex_cross.selected_parent, population.genes, population.children, population.elite_gene);
		/*-------------------------------------------------------*/

		/* エリート情報書き出し */		
		fpX.open("GA_adapt.csv", std::ios::app);
		population.output_elite(gg, total_eval_count, fpX);
		fpX.close();
		/* ５世代おきに集団情報書き出し */
		if(gg % 5 == 0){
			population.output_results(gg);
		}
	}
	/*-------------------------------------------------------*/
	population.output_results(GENERATION_SIZE );
}

/*//===========================================================
// ● 個体群の評価を実施評価
//=========================================================== */
void SingleObjOptimizer::eval_new_genes(vector<RGeneS>& rgenes){
	/* 定数群 */
	const int PARA_OMP_CALC_NUM = OptSettings::PARA_OMP_CALC_NUM;
	const int GENE_SIZE = OptSettings::GENE_SIZE;
	const int RGENE_FEATURES = OptSettings::get_total_oracle_features();
	/* このCallで評価する総個体数 */
	const int TOTAL_EVAL_NUM = rgenes.size();
	/* 並列分考慮したループ数 */
	const int LOOP_NUM1  = TOTAL_EVAL_NUM / PARA_OMP_CALC_NUM;
	const int LOOP_NUM2  = TOTAL_EVAL_NUM % PARA_OMP_CALC_NUM;

	/* 遺伝子形状 */
	double **gene_values = new double*[PARA_OMP_CALC_NUM];
	for(int i = 0 ; i < PARA_OMP_CALC_NUM ; i++){
		gene_values[i] = new double[GENE_SIZE];
	}
	/* 結果情報 */
	int* conv_flag = new int[PARA_OMP_CALC_NUM];
	double** oracle_results = new double*[PARA_OMP_CALC_NUM];
	for(int i = 0 ; i < PARA_OMP_CALC_NUM ; i++){
		oracle_results[i] = new double[RGENE_FEATURES];
	}
	
	/*-------------------------------------------------------*/
	/* ループその１・・・割り切れた分だけループを回す */
	for(int kkk = 0 ; kkk < LOOP_NUM1 ; kkk++){
		//cout << "Loop " << kkk << endl;
		for(int i = 0 ; i < PARA_OMP_CALC_NUM ; i++){
			const int pp = i + kkk*PARA_OMP_CALC_NUM;
			/* 遺伝子情報をセット */
			rgenes[pp].getGene(gene_values[i]);
			oracles[i].setGeneParas(gene_values[i]);
		}
		/* 評価を並列計算 */
#ifdef OMP_USING_OPT_EVAL
		#pragma omp parallel for
#endif
		for(int i = 0 ; i < PARA_OMP_CALC_NUM; i++){
			conv_flag[i] = oracles[i].run_oracle(oracle_results[i]);
		}
		/* 適応度情報を集団に渡す */
		for(int i = 0 ; i < PARA_OMP_CALC_NUM ; i++){
			const int pp = i + kkk*PARA_OMP_CALC_NUM;
			const int conv_type = (conv_flag[i] ? 1 : 0);
			rgenes[pp].setAdapts(oracle_results[i], conv_type);
			oracles[i].reset_oracle_data();
		}
	}
	/*-------------------------------------------------------*/
	/* ループその２・・・余った分 */
	if(LOOP_NUM2 != 0){
		for(int i = 0 ; i < LOOP_NUM2 ; i++){
			//cout << "Loop2 " << i << endl;
			const int pp = i + LOOP_NUM1*PARA_OMP_CALC_NUM;
			/* 遺伝子情報をセット */
			rgenes[pp].getGene(gene_values[i]);
			oracles[i].setGeneParas(gene_values[i]);
		}
		/* 評価を並列計算 */
#ifdef OMP_USING_OPT_EVAL
		#pragma omp parallel for
#endif
		for(int i = 0 ; i < LOOP_NUM2 ; i++){
			conv_flag[i] = oracles[i].run_oracle(oracle_results[i]);
		}
		/* 適応度情報を集団に渡す */
		for(int i = 0 ; i < LOOP_NUM2 ; i++){
			const int pp = i + LOOP_NUM1*PARA_OMP_CALC_NUM;
			const int conv_type = (conv_flag[i] ? 1 : 0);
			rgenes[pp].setAdapts(oracle_results[i], conv_type);
			oracles[i].reset_oracle_data();
		}
	}

	/* 評価回数カウント */
	this->total_eval_count += LOOP_NUM1 + LOOP_NUM2;
	
	for(int i = 0 ; i < PARA_OMP_CALC_NUM ; i++){
		delete[] gene_values[i];
	}
	delete[] gene_values;
	delete[] conv_flag;
	for(int i = 0 ; i < PARA_OMP_CALC_NUM ; i++){
		delete[] oracle_results[i];
	}
	delete[] oracle_results;
}


