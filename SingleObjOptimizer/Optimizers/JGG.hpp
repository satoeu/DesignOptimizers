#ifndef DEF_NSGA_III_DEF
#define DEF_NSGA_III_DEF

#include <BaseFuncs/Mt.hpp>
#include "../RGeneS.hpp"
#include "../RPopulationS.hpp"
#include "RootOptimizerS.hpp"


/*
//=============================================================
// ■ 世代交代モデルJGG
//=============================================================
//	
//============================================================= */
class JGG : public RootOptimizerS{
protected:
	/*--------------------------------------------------------------------*/
public:
	JGG(int p_num);												                    /* コンストラクタ */
	~JGG(){;};																		/* デストラクタ */
	/* 生存選択 */
    void selection(const vector<int>& selected_parent, vector<RGeneS>& ParentPop, vector<RGeneS>& Children, RGeneS& elite_gene);
};





#endif
