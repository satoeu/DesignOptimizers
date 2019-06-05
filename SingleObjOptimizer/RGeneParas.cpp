#include "RGeneS.hpp"
#include <BaseFuncs/Mt.hpp>

/*//===========================================================
  // ● 遺伝子座の値の決定
  //=========================================================== */
void RGeneS::makeGene(){
	for(int kkk = 0 ; kkk < 1500 ; kkk++){
		/* ランダムに遺伝子セット */
		for(int i=0 ; i< size ; i++){
			gene[i] = CommonLibs::Mt::unif_rand(-6.0, 6.0);
		}
		/* 遺伝子補正 */
		this->adjust();

		/* 形状生成 */
		//make_shape();
		
		bool ok = check_simple_violate_init();
		if(ok){
			break;
		}
	}
}

/*//===========================================================
  // ● 簡易的に見れる制約をチェック(初期個体)
  //=========================================================== */
bool RGeneS::check_simple_violate_init(){
	return(true);
}
/*//===========================================================
// ● 簡易的に見れる制約をチェック(子個体)
//=========================================================== */
bool RGeneS::check_simple_violate_child(){
	bool bl = check_simple_violate_init();
	return(bl);
}

/*//===========================================================
  // ● 遺伝子補正
  //=========================================================== */
void RGeneS::adjust(){
	/*
	for(int i=0 ; i< size ; i++){
		if(gene[i] < 0){
			gene[i] = 0;
		}else if(gene[i] >= 1.0){
			gene[i] = 1.0;
		}
	}*/
}

