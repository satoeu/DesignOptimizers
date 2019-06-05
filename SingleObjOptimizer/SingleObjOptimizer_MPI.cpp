#include "SingleObjOptimizer.hpp"
#include <fstream>

#ifdef MPI_RUN_OPT_MODE_SETTINGS

/*
//=============================================================
// ■ SingleObjOptimizer
//=============================================================
//	最適化実行主体クラス
//============================================================= */

/*//===========================================================
// ● 最適化実行
//=========================================================== */
void SingleObjOptimizer::run_optimization_MPI(){
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

	int mpi_rank, mpi_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

	fstream fpX;
	if(mpi_rank == 0){
		fstream fpX;
		fpX.open("GA_adapt.csv", std::ios::out);
		fpX << "generation, num_eval, adapt, ";
		for(int k=0 ; k < OptSettings::NUM_FEATURE ; k++){
			fpX << OptSettings::FEATURES_NAME[k] << ", ";
		}
		fpX << "residual" << endl;
		fpX.close();
	}

	/* OpenMP並列数セット */
	omp_set_num_threads(1);

	/* 評価数の統一（個体はルートノードでしか作らないから） */
	int eval_size = population.genes.size();
	int* eval_size_all = new int[PARA_OMP_CALC_NUM];
	for(int i = 0 ; i < PARA_OMP_CALC_NUM ; i++){
		eval_size_all[i] = eval_size;
	}
	int* eval_size_one = new int[1];
	int ierr1 = MPI_Scatter(eval_size_all, 1, MPI_INT, eval_size_one, 1, MPI_INT, 0, MPI_COMM_WORLD);
	/* 初期個体群評価 */
	this->eval_new_genes_MPI(*eval_size_one, population.genes);

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
		/*-------------------------------------------------------*/
		/* REX */
		if(mpi_rank == 0){
			cout << "Generation " << gg << endl;
			rex_cross.create_children(population.genes, population.children);
		}
		/* 子個体評価 */
		eval_size = population.children.size();
		for(int i = 0 ; i < PARA_OMP_CALC_NUM ; i++){
			eval_size_all[i] = eval_size;
		}
		ierr1 = MPI_Scatter(eval_size_all, 1, MPI_INT, eval_size_one, 1, MPI_INT, 0, MPI_COMM_WORLD);
		this->eval_new_genes_MPI(*eval_size_one, population.children);
		/*-------------------------------------------------------*/
		if(mpi_rank == 0){
			/* JGGで個体群更新 */
			jgg.selection(rex_cross.selected_parent, population.genes, population.children, population.elite_gene);
		}
		/*-------------------------------------------------------*/

		/* 簡易HV書き出し */
		if(mpi_rank == 0){
			/* エリート情報書き出し */		
			fpX.open("GA_adapt.csv", std::ios::app);
			population.output_elite(gg, total_eval_count, fpX);
			fpX.close();
			/* ５世代おきに集団情報書き出し */
			if(gg % 5 == 0){
				population.output_results(gg);
			}
		}
	}
	/*-------------------------------------------------------*/
	if(mpi_rank == 0){
		population.output_results(GENERATION_SIZE );
	}

	delete[] eval_size_all;
	delete[] eval_size_one;
}

/*//===========================================================
// ● 個体群を評価
//=========================================================== */
void SingleObjOptimizer::eval_new_genes_MPI(const int eval_size, vector<RGeneS>& rgenes){
	/* 定数群 */
	const int GENE_SIZE = OptSettings::GENE_SIZE;
	const int RGENE_FEATURES = OptSettings::get_total_oracle_features();
	/* 並列数（今の場合、MPI並列数） */
	const int PARA_OMP_CALC_NUM = OptSettings::PARA_OMP_CALC_NUM;

	/* このCallで評価する総個体数 */
	//const int TOTAL_EVAL_NUM = rgenes.size();
	/* 並列分考慮したループ数 */
	const int LOOP_NUM1  = eval_size / PARA_OMP_CALC_NUM;
	const int LOOP_NUM2  = eval_size % PARA_OMP_CALC_NUM;

	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	//cout << "SIZE " << mpi_rank << " >> "<< eval_size << ", " << LOOP_NUM1 << ", " << LOOP_NUM2 << endl;
	/* 遺伝子形状 */
	double *gene_values_all = new double[PARA_OMP_CALC_NUM*GENE_SIZE];
	double **gene_values = new double*[PARA_OMP_CALC_NUM];
	for(int i = 0 ; i < PARA_OMP_CALC_NUM ; i++){
		gene_values[i] = gene_values_all + i*GENE_SIZE;
	}
	double *gene_values_one = new double[GENE_SIZE];/*１遺伝子情報の受信用*/

	/* 結果フラグ */
	int* conv_flag_all = new int[PARA_OMP_CALC_NUM];
	int* conv_flag_one = new int[1];
	/* 結果情報 */
	double* oracle_results_all = new double[PARA_OMP_CALC_NUM*RGENE_FEATURES];
	double** oracle_results = new double*[PARA_OMP_CALC_NUM];
	for(int i = 0 ; i < PARA_OMP_CALC_NUM ; i++){
		oracle_results[i] = oracle_results_all + i*RGENE_FEATURES;
	}
	double* oracle_results_one = new double[RGENE_FEATURES];/*１遺伝子情報の受信用*/

	/*-------------------------------------------------------*/
	/* ループその１・・・割り切れた分だけループを回す */
	for(int kkk = 0 ; kkk < LOOP_NUM1 ; kkk++){
		/* ルートノードで、全遺伝子を１つの配列にまとめる */
		if(mpi_rank == 0){
			for(int i = 0 ; i < PARA_OMP_CALC_NUM ; i++){
				const int pp = i + kkk*PARA_OMP_CALC_NUM;
				rgenes[pp].getGene(gene_values[i]);
			}
		}
		/**/
		/* ルートノードから、全ノードへ遺伝子情報を１つずつ送信する */
		int ierr1 = MPI_Scatter(gene_values_all, GENE_SIZE, MPI_DOUBLE, gene_values_one, GENE_SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		/**/
		/* 遺伝子をOracleに代入 */
		oracles->setGeneParas(gene_values_one);

		/**/
		/* 評価を計算 */
		conv_flag_one[0] = oracles->run_oracle(oracle_results_one);
		/**/

		/* 結果をルートノードに集約する */
		int ierr2 = MPI_Gather(oracle_results_one, RGENE_FEATURES, MPI_DOUBLE, oracle_results_all, RGENE_FEATURES, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		/* 成功情報も */
		int ierr3 = MPI_Gather(conv_flag_one, 1, MPI_INT, conv_flag_all, 1, MPI_INT, 0, MPI_COMM_WORLD);

		/* 適応度情報を集団に渡す */
		if(mpi_rank == 0){
			for(int i = 0 ; i < PARA_OMP_CALC_NUM ; i++){
				const int pp = i + kkk*PARA_OMP_CALC_NUM;
				const int conv_type = (conv_flag_all[i] ? 1 : 0);
				rgenes[pp].setAdapts(oracle_results[i], conv_type);
				oracles[i].reset_oracle_data();
			}
		}
	}
	/*-------------------------------------------------------*/
	/* ループその２・・・余った分 */
	if(LOOP_NUM2 != 0){
		/* ルートノードで、全遺伝子を１つの配列にまとめる */
		if(mpi_rank == 0){
			for(int i = 0 ; i < LOOP_NUM2 ; i++){
				const int pp = i + LOOP_NUM1*PARA_OMP_CALC_NUM;
				rgenes[pp].getGene(gene_values[i]);
			}
		}
		/**/
		/* ルートノードから、全ノードへ遺伝子情報を１つずつ送信する(後ろは余るが無視すると) */
		int ierr1 = MPI_Scatter(gene_values_all, GENE_SIZE, MPI_DOUBLE, gene_values_one, GENE_SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		/**/
		/* 遺伝子をOracleに代入 */
		oracles->setGeneParas(gene_values_one);

		/**/
		/* 評価を計算 */
		if(mpi_rank < LOOP_NUM2){
			conv_flag_one[0] = oracles->run_oracle(oracle_results_one);
		}
		/**/

		/* 結果をルートノードに集約する */
		int ierr2 = MPI_Gather(oracle_results_one, RGENE_FEATURES, MPI_DOUBLE, oracle_results_all, RGENE_FEATURES, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		/* 成功情報も */
		int ierr3 = MPI_Gather(conv_flag_one, 1, MPI_INT, conv_flag_all, 1, MPI_INT, 0, MPI_COMM_WORLD);

		/* 適応度情報を集団に渡す */
		if(mpi_rank == 0){
			for(int i = 0 ; i < LOOP_NUM2 ; i++){
				const int pp = i + LOOP_NUM1*PARA_OMP_CALC_NUM;
				const int conv_type = (conv_flag_all[i] ? 1 : 0);
				rgenes[pp].setAdapts(oracle_results[i], conv_type);
				oracles[i].reset_oracle_data();
			}
		}
	}
	/* 評価回数カウント */
	this->total_eval_count += LOOP_NUM1 + LOOP_NUM2;

	delete[] gene_values;
	delete[] gene_values_all;
	delete[] gene_values_one;
	delete[] conv_flag_all;
	delete[] conv_flag_one;
	delete[] oracle_results;
	delete[] oracle_results_all;
	delete[] oracle_results_one;
}


#endif
