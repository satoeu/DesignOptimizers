#ifndef DEF_BASEDEF00_HEADER_DEF
#define DEF_BASEDEF00_HEADER_DEF

/* 整数型の定義～大規模にする際はlongにする */
using myint = int;

/* その他、世界定数 */
namespace CommonDef{
	constexpr double PI = 3.14159265358979;			/* 円周率 */
	constexpr double MYU0 = 4.0*PI*1.0e-7;			/* 真空透磁率 */
	constexpr double VNYU0 = 1.0/MYU0;				/* 真空磁気抵抗率 */
	constexpr double NORMB_EPS = 1.0e-10;			/* 微小ゼロとみなす磁束密度ノルム */
};

/* プログラム中のprintレベル制御 */
#define PRINTOUT_LEVEL1
#define PRINTOUT_LEVEL2

/* 実行環境がWindowsかどうか */
//#define IS_WINDOWS_SISTEM


#endif
