#ifndef DEF_MULTI_POPULATION_00
#define DEF_MULTI_POPULATION_00

#define _CRT_SECURE_NO_WARNINGS


#include <BaseDefines.hpp>
#include <vector>
#include <map>
#include "01_OptSetting.hpp"
#include "RGeneS.hpp"

/*
//=============================================================
// ■ RPopulation
//=============================================================
//	多数目的最適化の個体集団管理クラス
//============================================================= */
class RPopulationS{
public:
	int popu_size;												/* 集団のサイズ */
    vector<RGeneS> genes;										/* 遺伝子集合 */
    vector<RGeneS> children;									/* 子個体管理配列 */
    RGeneS elite_gene;				                            /* エリートアーカイブ */
/*===============================================================*/
	static void quicksort(int *GeneValues, const vector<RGeneS>& totalGenes, int up, int lo);		/* 適用度の順にクイックソート */
/*===============================================================*/
	RPopulationS(){return;};		                   							    /* コンストラクタ */
	RPopulationS(int popu_size00);									            /* コンストラクタ */
	~RPopulationS();
	void initialize(int popu_size00);									        /* コンストラクタ */
	int getPopuSize() const{return(genes.size());};
	int getChildSize() const{return(children.size());};
	void getParas(int i, double *x);											/* 個体ｉ番目の遺伝子データを配列ｘにコピー */
	void getChildParas(int i, double *x);										/* 子ｉ番目の遺伝子データを配列ｘにコピー */
	void calcAdapts(int i, double* oracle_results, int conv_type);				/* i番目の個体の適応度計算 */
	void calcChildAdapts(int i, double* oracle_results, int conv_type);			/* 子遺伝子i番目の個体の適応度計算 */
	/*--------------------------------------------------------------------*/    
    void output_elite(int gg, int num_e, fstream& fp);		                      			/* エリート情報書き出し */
	void output_results(int gen);
};

#endif
