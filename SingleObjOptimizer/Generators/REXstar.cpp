
#include "REXstar.hpp"

/*
//=============================================================
// ■ REXstar
//=============================================================
//	REXstar交叉による子個体生成（ただし範囲拡張などはない簡単化ver）
//============================================================= */


/*//===========================================================
// ● 子個体生成
//=========================================================== */
void REXstar::create_children(const vector<RGeneS>& populations, vector<RGeneS>& children){
	/* 一度親候補をまとめる */
	const int popu_size = populations.size();
	int* parents = new int[parent_num];
	/* 親をランダム抽出 */
	selectNewParents(parents, popu_size);

	/* 交叉開始！ */
	children.clear();
	crossREX(children, populations, parents);

	/* 選択された親を保存 */
	selected_parent.clear();
	for(int i = 0 ; i < parent_num ; i++){
		selected_parent.push_back(parents[i]);
	}

	delete[] parents;
}

/*//===========================================================
// ● 親個体の選択
//=========================================================== */
void REXstar::selectNewParents(int* parents, int total_size){
	/* 交叉に使う親選択 */
	for(int i = 0 ; i < parent_num ; i++){
		bool bl;
		int temp_i;
		/* 同じ親を選ばないようにループ */
		for(int kkk=0; kkk < 1000; kkk++){
			bl=false;
			temp_i = CommonLibs::Mt::mrand(total_size);
			for(int j = 0 ; j < i ; j++){
				bl |= (parents[j] == temp_i);
			}
			if(!bl) break;
		}
		/* 親に設定 */
		parents[i] = temp_i;
	}
}

/*//===========================================================
// ● REX交叉
//=========================================================== */
void REXstar::crossREX(vector<RGeneS>& children, const vector<RGeneS>& populations, int* parents){
	const int GENE_SIZE = OptSettings::GENE_SIZE;
	double *grav = new double[GENE_SIZE];
	/* 重心計算 */
	double tmp = 0;
	for(int i=0 ; i< GENE_SIZE ; i++){
		tmp = 0;
		for(int j=0 ; j< parent_num ; j++){
			const int par_pos = parents[j];
			tmp += populations[par_pos].getGene(i);
		}
		grav[i] = tmp / (1.0*parent_num);
	}

	/* 最良個体を選択 */
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

	/* REXアルゴリズム */
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
// ● 子個体生成メソッド部分
//=========================================================== */
void REXstar::child_generationREX(RGeneS& child, const vector<RGeneS>& populations, int* parents, double *par_grav, double *best_genes){
	const int GENE_SIZE = OptSettings::GENE_SIZE;
	double *rand_vec = new double[parent_num];
	double *child_gene = new double[GENE_SIZE];
	for(int xxx=0; xxx < 150 ; xxx++){
		/* 乱数生成 */
		for(int j = 0 ; j < parent_num ; j++){
			rand_vec[j] = CommonLibs::Mt::genrand_real1();
			rand_vec[j] -= 0.5;
			rand_vec[j] *= 2.0*sqrt(3.0/(parent_num));
		}
		double rand_1 = CommonLibs::Mt::genrand_real1();
		rand_1 *= rex_t;
		/* 子の遺伝子生成 */
		for(int j = 0 ; j < GENE_SIZE ; j++){
			/* REX部 */
			double temp_vec=0;
			for(int k = 0 ; k < parent_num ; k++){
				const int par_pos = parents[k];
				temp_vec += rand_vec[k] * (populations[par_pos].getGene(j) - par_grav[j]);
			}
			child_gene[j] = par_grav[j] + temp_vec;
			/* 降下方向 */
			child_gene[j] += rand_1 * (best_genes[j]-par_grav[j]);
		}
		/* 子を代入 */
		child.setGene(child_gene);
		bool ok = child.check_simple_violate_child();
		if(ok){
			break;
		}
	}
	delete[] rand_vec;
	delete[] child_gene;
}

