#include <fstream>
#include <BaseFuncs/Mt.hpp>
#include "RPopulationS.hpp"

/*
//=============================================================
// ■ RPopulationS
//=============================================================
//	多数目的最適化の個体集団管理クラス
//============================================================= */

/*//===========================================================
// ● コンストラクタ
//=========================================================== */
RPopulationS::RPopulationS(int popu_size00){
	this->initialize(popu_size00);
}

void RPopulationS::initialize(int popu_size00){
    popu_size = popu_size00;
    genes.resize(popu_size);
    for(int i = 0 ; i < popu_size ; i++){
        genes[i].makeGene();
    }
	children.clear();
	elite_gene.setAdaptNormal(1.0e+12);
}
/*//===========================================================
// ● デストラクタ
//=========================================================== */
RPopulationS::~RPopulationS(){
    genes.clear();
    children.clear();
}


/*//===========================================================
  // ● 個体ｉ番目の遺伝子データを配列ｘにコピー
  //=========================================================== */
void RPopulationS::getParas(int i, double *x){
	genes[i].getGene(x);
}

/*//===========================================================
  // ● 子ｉ番目の遺伝子データを配列ｘにコピー
  //=========================================================== */
void RPopulationS::getChildParas(int i, double *x){
	children[i].getGene(x);
}

/*//===========================================================
  // ● i番目の個体の適応度計算 
  //=========================================================== */
void RPopulationS::calcAdapts(int i, double* oracle_results, int conv_type){
	genes[i].setAdapts(oracle_results, conv_type);
}

/*//===========================================================
  // ● i子遺伝子i番目の個体の適応度計算
  //=========================================================== */
void RPopulationS::calcChildAdapts(int i, double* oracle_results, int conv_type){
	children[i].setAdapts(oracle_results, conv_type);
}

/*//===========================================================
// ● 評価順にソート
//=========================================================== */
void RPopulationS::quicksort(int *GeneValues, const vector<RGeneS>& totalGenes, int up, int lo){
	int mid = (lo + up) / 2;
	double bound = totalGenes[GeneValues[mid]].getAdapts();
	int l = lo, u = up;
	do{
		double lV = totalGenes[GeneValues[l]].getAdapts();
		double uV = totalGenes[GeneValues[u]].getAdapts();
		while(lV < bound){
			l++;
			lV = totalGenes[GeneValues[l]].getAdapts();
		}
		while(uV > bound){
			u--;
			uV = totalGenes[GeneValues[u]].getAdapts();
		}
		if(l <= u){
			int temp = GeneValues[l];
			GeneValues[l] = GeneValues[u];
			GeneValues[u] = temp;
			l++;
			u--;
		}
	}while(l < u);
	if(lo < u) quicksort(GeneValues, totalGenes, u, lo);
	if(l < up) quicksort(GeneValues, totalGenes, up, l);
}



/*//===========================================================
// ● エリート情報書き出し
//=========================================================== */
void RPopulationS::output_elite(int gg, int num_e, fstream& fp){
	double ad = elite_gene.getAdapts();
	fp << gg << ", " << num_e << ", " << ad << ", ";
	for(int k=0 ; k < OptSettings::NUM_FEATURE ; k++){
		double temp = elite_gene.getFeatureVal(k);
		fp << temp << ", " ;
	}
	ad = elite_gene.getResidual();
	fp << ad << endl;
}



/*//===========================================================
  // ● 結果書き出し
  //=========================================================== */
void RPopulationS::output_results(int gen){
	string str1 = "./Data/Results";
	string str2 = ".csv";
	string str_ID = std::to_string(gen);
	string filename =str1 + str_ID + str2;
//	system("rm ./Data/*");
	fstream fp1(filename, std::ios::out);

	const int current_pop_size = genes.size();	
	fp1 << "ID,adapt,";
	for(int k=0 ; k < OptSettings::NUM_FEATURE ; k++){
		fp1 << OptSettings::FEATURES_NAME[k] << ", ";
	}
	fp1 << "residual"<< endl;;

	for(int i=0 ; i< current_pop_size ; i++){
		/* 目的関数値を出力 */
		fp1 << i << "," << genes[i].getAdapts() << ",";
		for(int k=0 ; k < OptSettings::NUM_FEATURE ; k++){
			fp1 << scientific << genes[i].getFeatureVal(k) << ", ";
		}
		fp1 << scientific <<genes[i].getResidual() << endl;
		/* 各個体の詳細情報書き出し */
		/*string str11 = "./Data/Gene/ResultGene";
		string str12 = ".csv";
		string str_ID2 = std::to_string(i);
		string filename2 =str11 + str_ID2 + str12;
		genes[i].PrintGene(filename2);*/
		/*str11 = "./Data/Gene/ResultGene";
		str12 = ".json";
		str_ID2 = std::to_string(i);
		filename2 =str11 + str_ID2 + str12;
		genes[i].PrintGeneJSON(filename2);*/
	}
	elite_gene.PrintGene("./Data/Gene/ResultsElite.csv");
	elite_gene.PrintGeneJSON("./Data/Gene/ResultsElite.json");
	fp1.close();
}


