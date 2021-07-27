#ifndef SMMAP_H
#define SMMAP_H
#include <gmapping/grid/map.h>
#include <gmapping/grid/harray2d.h>
#include <gmapping/utils/point.h>
#define SIGHT_INC 1

namespace GMapping {
//将击中的点位置叠加,并且求均值
struct PointAccumulator{
	typedef point<float> FloatPoint;
	/* before 
	PointAccumulator(int i=-1): acc(0,0), n(0), visits(0){assert(i==-1);}
	*/
	/*after begin*/
	PointAccumulator(): acc(0,0), n(0), visits(0){}
	PointAccumulator(int i): acc(0,0), n(0), visits(0){assert(i==-1);}
	/*after end*/
    inline void update(bool value, const Point& p=Point(0,0));
	inline Point mean() const {return 1./n*Point(acc.x, acc.y);}//这个均值是不是就是机器人的位姿
	inline operator double() const { return visits?(double)n*SIGHT_INC/(double)visits:-1; }
	inline void add(const PointAccumulator& p) {acc=acc+p.acc; n+=p.n; visits+=p.visits; }//p表示激光束经过栅格坐标吗？？？？？？？？？？？？？？？？？坐标叠加
	static const PointAccumulator& Unknown();
	static PointAccumulator* unknown_ptr;
	FloatPoint acc;// //坐标累计值
	int n, visits;//n激光hit的次数,visits 被访问的次数
	inline double entropy() const;
};

void PointAccumulator::update(bool value, const Point& p){
	if (value) {////ture:更新    false:访问
		acc.x+= static_cast<float>(p.x);
		acc.y+= static_cast<float>(p.y); 
		n++; 
		visits+=SIGHT_INC;
	} else
		visits++;
}

double PointAccumulator::entropy() const{//求熵(entropy)的函数,通过记录累积的次数和访问的次数计算熵值
	if (!visits)
		return -log(.5);
	if (n==visits || n==0)
		return 0;
	double x=(double)n*SIGHT_INC/(double)visits;//通过记录累积的次数和访问的次数计算熵值
	return -( x*log(x)+ (1-x)*log(1-x) );
}


typedef Map<PointAccumulator,HierarchicalArray2D<PointAccumulator> > ScanMatcherMap;

};

#endif 
