
#include "REXstar.hpp"

/*
//=============================================================
// �� REXstar
//=============================================================
//	REXstar�����ɂ��q�̐����i�������͈͊g���Ȃǂ͂Ȃ��ȒP��ver�j
//============================================================= */


/*//===========================================================
// �� �q�̐���
//=========================================================== */
void REXstar::create_children(const vector<RGeneS>& populations, vector<RGeneS>& children){
	/* ��x�e�����܂Ƃ߂� */
	const int popu_size = populations.size();
	int* parents = new int[parent_num];
	/* �e�������_�����o */
	selectNewParents(parents, popu_size);

	/* �����J�n�I */
	children.clear();
	crossREX(children, populations, parents);

	/* �I�����ꂽ�e��ۑ� */
	selected_parent.clear();
	for(int i = 0 ; i < parent_num ; i++){
		selected_parent.push_back(parents[i]);
	}

	delete[] parents;
}

/*//===========================================================
// �� �e�̂̑I��
//=========================================================== */
void REXstar::selectNewParents(int* parents, int total_size){
	/* �����Ɏg���e�I�� */
	for(int i = 0 ; i < parent_num ; i++){
		bool bl;
		int temp_i;
		/* �����e��I�΂Ȃ��悤�Ƀ��[�v */
		for(int kkk=0; kkk < 1000; kkk++){
			bl=false;
			temp_i = CommonLibs::Mt::mrand(total_size);
			for(int j = 0 ; j < i ; j++){
				bl |= (parents[j] == temp_i);
			}
			if(!bl) break;
		}
		/* �e�ɐݒ� */
		parents[i] = temp_i;
	}
}

/*//===========================================================
// �� REX����
//=========================================================== */
void REXstar::crossREX(vector<RGeneS>& children, const vector<RGeneS>& populations, int* parents){
	const int GENE_SIZE = OptSettings::GENE_SIZE;
	double *grav = new double[GENE_SIZE];
	/* �d�S�v�Z */
	double tmp = 0;
	for(int i=0 ; i< GENE_SIZE ; i++){
		tmp = 0;
		for(int j=0 ; j< parent_num ; j++){
			const int par_pos = parents[j];
			tmp += populations[par_pos].getGene(i);
		}
		grav[i] = tmp / (1.0*parent_num);
	}

	/* �ŗǌ̂�I�� */
	double *best_genes = new double[GENE_SIZE];
	double best = 1.0e+12;
	int best_gp = 0;
	for(int i = 0 ; i < parent_num ; i++){
		const int par_pos = parents[i];
		double xxx = populations[par_pos].getAdapts();
		if(xxx < best){
			best = xxx;
			best_gp = i;
		}
	}
	for(int i = 0 ; i < GENE_SIZE ; i++){
		const int par_pos = parents[best_gp];
		best_genes[i] = populations[par_pos].getGene(i);
	}

	/* REX�A���S���Y�� */
	children.clear();
	for(int ccc = 0 ; ccc < child_num ; ccc++){
		RGeneS temp_child;
		child_generationREX(temp_child, populations, parents, grav, best_genes);
		children.push_back( std::move(temp_child) );
	}
	delete[] grav;
	delete[] best_genes;
}

/*//===========================================================
// �� �q�̐������\�b�h����
//=========================================================== */
void REXstar::child_generationREX(RGeneS& child, const vector<RGeneS>& populations, int* parents, double *par_grav, double *best_genes){
	const int GENE_SIZE = OptSettings::GENE_SIZE;
	double *rand_vec = new double[parent_num];
	double *child_gene = new double[GENE_SIZE];
	for(int xxx=0; xxx < 150 ; xxx++){
		/* �������� */
		for(int j = 0 ; j < parent_num ; j++){
			rand_vec[j] = CommonLibs::Mt::genrand_real1();
			rand_vec[j] -= 0.5;
			rand_vec[j] *= 2.0*sqrt(3.0/(parent_num));
		}
		double rand_1 = CommonLibs::Mt::genrand_real1();
		rand_1 *= rex_t;
		/* �q�̈�`�q���� */
		for(int j = 0 ; j < GENE_SIZE ; j++){
			/* REX�� */
			double temp_vec=0;
			for(int k = 0 ; k < parent_num ; k++){
				const int par_pos = parents[k];
				temp_vec += rand_vec[k] * (populations[par_pos].getGene(j) - par_grav[j]);
			}
			child_gene[j] = par_grav[j] + temp_vec;
			/* �~������ */
			child_gene[j] += rand_1 * (best_genes[j]-par_grav[j]);
		}
		/* �q���� */
		child.setGene(child_gene);
		bool ok = child.check_simple_violate_child();
		if(ok){
			break;
		}
	}
	delete[] rand_vec;
	delete[] child_gene;
}

