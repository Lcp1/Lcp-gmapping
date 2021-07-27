#ifndef SCANMATCHER_H
#define SCANMATCHER_H

#include "gmapping/scanmatcher/icp.h"
#include "gmapping/scanmatcher/smmap.h"
#include <gmapping/utils/macro_params.h>
#include <gmapping/utils/stat.h>
#include <iostream>
#include <gmapping/utils/gvalues.h>
#define LASER_MAXBEAMS 2048

namespace GMapping {

class ScanMatcher{
	public:
		typedef Covariance3 CovarianceMatrix;
		
		ScanMatcher();
		~ScanMatcher();
		double icpOptimize(OrientedPoint& pnew, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const;
		double optimize(OrientedPoint& pnew, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const;
		double optimize(OrientedPoint& mean, CovarianceMatrix& cov, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const;
		
		double   registerScan(ScanMatcherMap& map, const OrientedPoint& p, const double* readings);
		void setLaserParameters
			(unsigned int beams, double* angles, const OrientedPoint& lpose);
		void setMatchingParameters
			(double urange, double range, double sigma, int kernsize, double lopt, double aopt, int iterations, double likelihoodSigma=1, unsigned int likelihoodSkip=0 );
		void invalidateActiveArea();
		void computeActiveArea(ScanMatcherMap& map, const OrientedPoint& p, const double* readings);

		inline double icpStep(OrientedPoint & pret, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const;
		inline double score(const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const;
		inline unsigned int likelihoodAndScore(double& s, double& l, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const;
		double likelihood(double& lmax, OrientedPoint& mean, CovarianceMatrix& cov, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings);
		double likelihood(double& _lmax, OrientedPoint& _mean, CovarianceMatrix& _cov, const ScanMatcherMap& map, const OrientedPoint& p, Gaussian3& odometry, const double* readings, double gain=180.);
		inline const double* laserAngles() const { return m_laserAngles; }
		inline unsigned int laserBeams() const { return m_laserBeams; }
		
		static const double nullLikelihood;
	protected:
		//state of the matcher
		bool m_activeAreaComputed;
		
		/**laser parameters*/
		unsigned int m_laserBeams;
		double       m_laserAngles[LASER_MAXBEAMS];
		//OrientedPoint m_laserPose;
		PARAM_SET_GET(OrientedPoint, laserPose, protected, public, public)
		PARAM_SET_GET(double, laserMaxRange, protected, public, public)
		/**scan_matcher parameters*/
		PARAM_SET_GET(double, usableRange, protected, public, public)
		PARAM_SET_GET(double, gaussianSigma, protected, public, public)
		PARAM_SET_GET(double, likelihoodSigma, protected, public, public)
		PARAM_SET_GET(int,    kernelSize, protected, public, public)
		PARAM_SET_GET(double, optAngularDelta, protected, public, public)
		PARAM_SET_GET(double, optLinearDelta, protected, public, public)
		PARAM_SET_GET(unsigned int, optRecursiveIterations, protected, public, public)
		PARAM_SET_GET(unsigned int, likelihoodSkip, protected, public, public)
		PARAM_SET_GET(double, llsamplerange, protected, public, public)
		PARAM_SET_GET(double, llsamplestep, protected, public, public)
		PARAM_SET_GET(double, lasamplerange, protected, public, public)
		PARAM_SET_GET(double, lasamplestep, protected, public, public)
		PARAM_SET_GET(bool, generateMap, protected, public, public)
		PARAM_SET_GET(double, enlargeStep, protected, public, public)
		PARAM_SET_GET(double, fullnessThreshold, protected, public, public)
		PARAM_SET_GET(double, angularOdometryReliability, protected, public, public)
		PARAM_SET_GET(double, linearOdometryReliability, protected, public, public)
		PARAM_SET_GET(double, freeCellRatio, protected, public, public)
		PARAM_SET_GET(unsigned int, initialBeamsSkip, protected, public, public)
		// allocate this large array only once
		IntPoint* m_linePoints;
};

inline double ScanMatcher::icpStep(OrientedPoint & pret, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const{
	const double * angle = m_laserAngles+m_initialBeamsSkip;
	OrientedPoint lp = p;
	lp.x += cos(p.theta)*m_laserPose.x-sin(p.theta)*m_laserPose.y;
	lp.y += sin(p.theta)*m_laserPose.x+cos(p.theta)*m_laserPose.y;
	lp.theta+=m_laserPose.theta;
	unsigned int skip=0;
	double freeDelta=map.getDelta()*m_freeCellRatio;
	std::list<PointPair> pairs;
	
	for (const double* r=readings+m_initialBeamsSkip; r<readings+m_laserBeams; r++, angle++){
		skip++;
		skip=skip>m_likelihoodSkip?0:skip;
		if (*r>m_usableRange||*r==0.0) continue;
		if (skip) continue;
		Point phit=lp;
		phit.x+=*r*cos(lp.theta+*angle);
		phit.y+=*r*sin(lp.theta+*angle);
		IntPoint iphit=map.world2map(phit);
		Point pfree=lp;
		pfree.x+=(*r-map.getDelta()*freeDelta)*cos(lp.theta+*angle);
		pfree.y+=(*r-map.getDelta()*freeDelta)*sin(lp.theta+*angle);
 		pfree=pfree-phit;
		IntPoint ipfree=map.world2map(pfree);
		bool found=false;
		Point bestMu(0.,0.);
		Point bestCell(0.,0.);
		for (int xx=-m_kernelSize; xx<=m_kernelSize; xx++)
		for (int yy=-m_kernelSize; yy<=m_kernelSize; yy++){
			IntPoint pr=iphit+IntPoint(xx,yy);
			IntPoint pf=pr+ipfree;
			//AccessibilityState s=map.storage().cellState(pr);
			//if (s&Inside && s&Allocated){
				const PointAccumulator& cell=map.cell(pr);
				const PointAccumulator& fcell=map.cell(pf);
				if (((double)cell )> m_fullnessThreshold && ((double)fcell )<m_fullnessThreshold){
					Point mu=phit-cell.mean();
					if (!found){
						bestMu=mu;
						bestCell=cell.mean();
						found=true;
					}else
						if((mu*mu)<(bestMu*bestMu)){
							bestMu=mu;
							bestCell=cell.mean();
						} 
						
				}
			//}
		}
		if (found){
			pairs.push_back(std::make_pair(phit, bestCell));
			//std::cerr << "(" << phit.x-bestCell.x << "," << phit.y-bestCell.y << ") ";
		}
		//std::cerr << std::endl;
	}
	
	OrientedPoint result(0,0,0);
	//double icpError=icpNonlinearStep(result,pairs);
	std::cerr << "result(" << pairs.size() << ")=" << result.x << " " << result.y << " " << result.theta << std::endl;
	pret.x=p.x+result.x;
	pret.y=p.y+result.y;
	pret.theta=p.theta+result.theta;
	pret.theta=atan2(sin(pret.theta), cos(pret.theta));
	return score(map, p, readings);
}
/*
* 计算每个激光点，在搜索框下面找到与初始点欧式距离最近的点，以最近距离计算分布。
* p当前位姿 
* readings点云数据
* map 地图
*/ 
inline double ScanMatcher::score(const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const{
	double s=0;
	const double * angle=m_laserAngles+m_initialBeamsSkip;//应该是点云的角度
	OrientedPoint lp=p;     //获取机器人位姿
	//先获取雷达在世界的坐标系
	lp.x+=cos(p.theta)*m_laserPose.x-sin(p.theta)*m_laserPose.y;
	lp.y+=sin(p.theta)*m_laserPose.x+cos(p.theta)*m_laserPose.y;
	//根据雷达坐标将点云变换到世界坐标系
	lp.theta+=m_laserPose.theta;
	unsigned int skip=0;
    /*
     * map.getDelta表示地图分辨率 m_freeCellRatio = sqrt(2)
     * 意思是如果激光hit了某个点 那么沿着激光方向的freeDelta距离的地方要是空闲才可以
    */
	double freeDelta=map.getDelta()*m_freeCellRatio;
	//枚举所有激光束,r为激光束的长度
	for (const double* r=readings+m_initialBeamsSkip; r<readings+m_laserBeams; r++, angle++){
		//设定跳过激光束
		skip++;
		skip=skip>m_likelihoodSkip?0:skip;//达到似然域跳过数量????????,就重置skip
		if (skip||*r>m_usableRange||*r==0.0) continue;//
		Point phit=lp;
		//计算出hit的点在世界的坐标
		phit.x+=*r*cos(lp.theta+*angle);
		phit.y+=*r*sin(lp.theta+*angle);
		//激光点世界坐标转到地图坐标
		IntPoint iphit=map.world2map(phit);
		Point pfree=lp;//激光束方向的空闲终点
		pfree.x+=(*r-map.getDelta()*freeDelta)*cos(lp.theta+*angle);//将雷达坐标加上hit前的空闲激光束坐标（相对激光雷达坐标）=空闲点的世界坐标
		pfree.y+=(*r-map.getDelta()*freeDelta)*sin(lp.theta+*angle);
 		pfree=pfree-phit;//空闲点终点到激光hit点的向量
		IntPoint ipfree=map.world2map(pfree);//转成地图坐标
		bool found=false;
		Point bestMu(0.,0.);
        /*在kernelSize大小的窗口中搜索出最优最可能被这个激光束击中的点 这个kernelSize在大小为1*/
		//m_kernelSize对应概率机器人p130的第七行m
		for (int xx=-m_kernelSize; xx<=m_kernelSize; xx++)
		for (int yy=-m_kernelSize; yy<=m_kernelSize; yy++){
			IntPoint pr=iphit+IntPoint(xx,yy);
			IntPoint pf=pr+ipfree;
			//AccessibilityState s=map.storage().cellState(pr);
			//if (s&Inside && s&Allocated){
				//根据坐标对应栅格,//PointAccumulator将击中的点位置叠加,并且求均值
				const PointAccumulator& cell=map.cell(pr);//记录击中的坐标值，直接引用cell，指向pr坐标的map.cell(pr)
				const PointAccumulator& fcell=map.cell(pf);//记录空闲坐标值
				if (((double)cell )> m_fullnessThreshold && ((double)fcell )<m_fullnessThreshold){//m_fullnessThreshold =0.1
				// 	求phit障碍点与均值的误差
					Point mu=phit-cell.mean();//mean 1./n*Point(acc.x, acc.y)其中cell的坐标计算如下，n为被击中的次数，acc.x和acc.y是击中的累加值：
					if (!found){//第一次时候将bestMu=mu
						bestMu=mu;
						found=true;
					}else
						bestMu=(mu*mu)<(bestMu*bestMu)?mu:bestMu;//选取距离最短的
				}
			//}
		}
		if (found)
		//高斯分布获取得分,对概率机器人p130的prob(dist,6hit)
			s+=exp(-1./m_gaussianSigma*bestMu*bestMu);
	}
	return s;
}

/**
* 计算某一个机器人的位置的似然，这个似然是在机器人位置的一个范围内( -mllsamplerange~mllsamplerange -mlasamplerange m_lasamplerange)的似然
* 并且通过这个范围内的似然分布，来求出机器人的位置的期望值和方差
* 这个函数的真实意图应该是计算机器人位姿p的似然值。
* 这个也是Cyrill Stachniss在论文里面提出的方法。
* 这个方法应该是和上面的optimize(_mean,_cov,...)方法配套使用的。
* 用optimize来计算位姿 然后用likelihoodAndScore函数来计算似然(或者说对应的位姿的权重)
* 可以认为这个函数的返回值就是对应的粒子的权重
* 这个函数的意思是：
* 假设机器人在p附近的
* (-mllsamplerange~mllsamplerange  于似然计算的平移采样距离
* -mlasamplerange m_lasamplerange)用于似然计算的角度采样距离
* 窗口是服从高斯分布的。
* 求机器人的真实位姿以及机器人的方差
@ s: 分布概率分数
@ l: 似然
@ p: 机器人位姿
@ readings:
*/

inline unsigned int ScanMatcher::likelihoodAndScore(double& s, double& l, const ScanMatcherMap& map, const OrientedPoint& p, const double* readings) const{
	using namespace std;
	l=0;//初始化参数
	s=0;
	const double * angle=m_laserAngles+m_initialBeamsSkip;//激光桢起始点云角度?
	OrientedPoint lp=p;//机器人位姿
	lp.x+=cos(p.theta)*m_laserPose.x-sin(p.theta)*m_laserPose.y;//计算出雷达位姿
	lp.y+=sin(p.theta)*m_laserPose.x+cos(p.theta)*m_laserPose.y;
	lp.theta+=m_laserPose.theta;//计算出雷达在世界的角度
	double noHit=nullLikelihood/(m_likelihoodSigma);//-.5 / 
	unsigned int skip=0;
	unsigned int c=0;
	double freeDelta=map.getDelta()*m_freeCellRatio;
	for (const double* r=readings+m_initialBeamsSkip; r<readings+m_laserBeams; r++, angle++){
		skip++;
		skip=skip>m_likelihoodSkip?0:skip;
		if (*r>m_usableRange) continue;
		if (skip) continue;
		Point phit=lp;
		phit.x+=*r*cos(lp.theta+*angle);
		phit.y+=*r*sin(lp.theta+*angle);
		IntPoint iphit=map.world2map(phit);//hit在世界的坐标
		Point pfree=lp;
		pfree.x+=(*r-freeDelta)*cos(lp.theta+*angle);
		pfree.y+=(*r-freeDelta)*sin(lp.theta+*angle);
		pfree=pfree-phit;//空闲的距离
		IntPoint ipfree=map.world2map(pfree);//空闲在地图的坐标   // 计算激光击中点前一个空白区域坐标
		bool found=false;
		Point bestMu(0.,0.);
		//在地图窗口下循环,计算最小距离值计算概率值最大值
		 // 初始激光点坐标基础下，在窗口下移动激光点坐标
		for (int xx=-m_kernelSize; xx<=m_kernelSize; xx++)
		for (int yy=-m_kernelSize; yy<=m_kernelSize; yy++){
			IntPoint pr=iphit+IntPoint(xx,yy);
			IntPoint pf=pr+ipfree;
			//AccessibilityState s=map.storage().cellState(pr);
			//if (s&Inside && s&Allocated){
				const PointAccumulator& cell=map.cell(pr);
				const PointAccumulator& fcell=map.cell(pf);
				if (((double)cell )>m_fullnessThreshold && ((double)fcell )<m_fullnessThreshold){
					// 计算欧式距离
					Point mu=phit-cell.mean();
					if (!found){
						bestMu=mu;
						found=true;
					}else
						bestMu=(mu*mu)<(bestMu*bestMu)?mu:bestMu;
				}
			//}	
		}
		if (found){
			 // 计算分布概率分数
			s+=exp(-1./m_gaussianSigma*bestMu*bestMu);//计算分值
			c++;
		}
		if (!skip){  /*计算似然 似然是计算所有的激光束 如果某一个激光束打中了空地 那也需要计算进去*/
			double f=(-1./m_likelihoodSigma)*(bestMu*bestMu);
			l+=(found)?f:noHit;
		}
	}
	return c;
}

};

#endif
