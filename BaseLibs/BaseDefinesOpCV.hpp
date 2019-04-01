#ifndef DEF_BASEDEF_HEADER_OPENCV_VER_DEF
#define DEF_BASEDEF_HEADER_OPENCV_VER_DEF

/* OpenCVを使うかどうか */
#define USING_CV_DRAW

/* OpenCVのlibの名前 */
#ifdef USING_CV_DRAW
#ifdef IS_WINDOWS_SISTEM
#define	OPENCV_LIV		"opencv_world345.lib"
#define	OPENCV_LIV_D	"opencv_world345d.lib"
#endif
#endif


#endif
