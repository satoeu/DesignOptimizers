
#include "JGG.hpp"

/*
//=============================================================
// ■ 世代交代モデルJGG
//=============================================================
//	
//============================================================= */

/*//===========================================================
  // ● コンストラクタ
  //===========================================================*/
JGG::JGG(int p_num){
	/* サイズセット */
	population_size = p_num;
}

/*//===========================================================
  // ● 生存選択
  //=========================================================== */
void JGG::selection(const vector<int>& selected_parent, vector<RGeneS>& ParentPop, vector<RGeneS>& Children, RGeneS& elite_gene){
	const int pop_size = ParentPop.size();
	const int child_num = Children.size();
	const int par_num = selected_parent.size();
	/* 交代する実際の数～少ないほうを優先 */
	const int sel_num_true = (par_num > child_num ? child_num : par_num);

	/* 評価順にソート */
	int* child_sort = new int[child_num];
	for(int i=0 ; i < child_num ; i++){
		child_sort[i] = i;
	}
	RPopulationS::quicksort(child_sort, Children, child_num-1, 0);

	/* エリート更新 */
	double best = elite_gene.getAdapts();
	for(int i=0 ; i < child_num ; i++){
		double ad = Children[i].getAdapts();		
		if(ad < best){
			best = ad;
			elite_gene = Children[i];
		}
	}
	/* 親と交代 */
	for(int i=0 ; i < sel_num_true ; i++){
		const int pos = selected_parent[i];
		ParentPop[pos] = Children[child_sort[i]];
	}
	delete[] child_sort;
}

