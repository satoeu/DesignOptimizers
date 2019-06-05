#ifndef DEF_REXSTAR00_CROSS_DEF
#define DEF_REXSTAR00_CROSS_DEF

#include <BaseFuncs/Mt.hpp>
#include "../RGeneS.hpp"
#include "RootCGeneratorS.hpp"


/*
//=============================================================
// �� REXstar
//=============================================================
//	REXstar�����ɂ��q�̐����i�������͈͊g���Ȃǂ͂Ȃ��ȒP��ver�j
//============================================================= */
class REXstar : public RootCGeneratorS{
protected:
	double rex_t;												/* �~�������X�e�b�v�T�C�Y */
	int parent_num;                                             /* �e�̂̐� */
/*--------------------------------------------------------------------*/
	void selectNewParents(int* parents, int total_size);																				/* �e�̂̃����_���I�� */
	void crossREX(vector<RGeneS>& children, const vector<RGeneS>& populations, int* parents);											/* REX���� */
	void child_generationREX(RGeneS& child, const vector<RGeneS>& populations, int* parents, double *par_grav, double *best_genes);		/* �������ԕ��� */
/*--------------------------------------------------------------------*/
public:
	vector<int> selected_parent;
	/* �R���X�g���N�^ */
	REXstar(int c_num, int p_num, double tt) : RootCGeneratorS(c_num){        
		child_num=c_num;
		parent_num = p_num;
		rex_t = tt;
	};
	~REXstar(){;};																	/* �f�X�g���N�^ */
																						/* �q�̐��� */
	void create_children(const vector<RGeneS>& populations, vector<RGeneS>& children);
};


#endif

