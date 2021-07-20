#include <string>
#include <deque>
#include <list>
#include <map>
#include <set>
#include <fstream>
#include <iomanip>
#include <gmapping/utils/stat.h>
#include "gmapping/gridfastslam/gridslamprocessor.h"

//#define MAP_CONSISTENCY_CHECK
//#define GENERATE_TRAJECTORIES

namespace GMapping {

const double m_distanceThresholdCheck = 20;
 
using namespace std;

  GridSlamProcessor::GridSlamProcessor(): m_infoStream(cout){
    
    period_ = 5.0;
    m_obsSigmaGain=1;
    m_resampleThreshold=0.5;
    m_minimumScore=0.;
  }
  
  GridSlamProcessor::GridSlamProcessor(const GridSlamProcessor& gsp) 
    :last_update_time_(0.0), m_particles(gsp.m_particles), m_infoStream(cout){
    
    period_ = 5.0;

    m_obsSigmaGain=gsp.m_obsSigmaGain;
    m_resampleThreshold=gsp.m_resampleThreshold;
    m_minimumScore=gsp.m_minimumScore;
    
    m_beams=gsp.m_beams;
    m_indexes=gsp.m_indexes;
    m_motionModel=gsp.m_motionModel;
    m_resampleThreshold=gsp.m_resampleThreshold;
    m_matcher=gsp.m_matcher;
    
    m_count=gsp.m_count;
    m_readingCount=gsp.m_readingCount;
    m_lastPartPose=gsp.m_lastPartPose;
    m_pose=gsp.m_pose;
    m_odoPose=gsp.m_odoPose;
    m_linearDistance=gsp.m_linearDistance;
    m_angularDistance=gsp.m_angularDistance;
    m_neff=gsp.m_neff;
	
    cerr << "FILTER COPY CONSTRUCTOR" << endl;
    cerr << "m_odoPose=" << m_odoPose.x << " " <<m_odoPose.y << " " << m_odoPose.theta << endl;
    cerr << "m_lastPartPose=" << m_lastPartPose.x << " " <<m_lastPartPose.y << " " << m_lastPartPose.theta << endl;
    cerr << "m_linearDistance=" << m_linearDistance << endl;
    cerr << "m_angularDistance=" << m_linearDistance << endl;
    
		
    m_xmin=gsp.m_xmin;
    m_ymin=gsp.m_ymin;
    m_xmax=gsp.m_xmax;
    m_ymax=gsp.m_ymax;
    m_delta=gsp.m_delta;
    
    m_regScore=gsp.m_regScore;
    m_critScore=gsp.m_critScore;
    m_maxMove=gsp.m_maxMove;
    
    m_linearThresholdDistance=gsp.m_linearThresholdDistance;
    m_angularThresholdDistance=gsp.m_angularThresholdDistance;
    m_obsSigmaGain=gsp.m_obsSigmaGain;
    
#ifdef MAP_CONSISTENCY_CHECK
    cerr << __PRETTY_FUNCTION__ <<  ": trajectories copy.... ";
#endif
    TNodeVector v=gsp.getTrajectories();
    for (unsigned int i=0; i<v.size(); i++){
		m_particles[i].node=v[i];
    }
#ifdef MAP_CONSISTENCY_CHECK
    cerr <<  "end" << endl;
#endif


    cerr  << "Tree: normalizing, resetting and propagating weights within copy construction/cloneing ..." ;
    updateTreeWeights(false);
    cerr  << ".done!" <<endl;
  }
  
  GridSlamProcessor::GridSlamProcessor(std::ostream& infoS): m_infoStream(infoS){
    period_ = 5.0;
    m_obsSigmaGain=1;
    m_resampleThreshold=0.5;
    m_minimumScore=0.;
	
  }

  GridSlamProcessor* GridSlamProcessor::clone() const {
# ifdef MAP_CONSISTENCY_CHECK
    cerr << __PRETTY_FUNCTION__ << ": performing preclone_fit_test" << endl;
    typedef std::map<autoptr< Array2D<PointAccumulator> >::reference* const, int> PointerMap;
    PointerMap pmap;
	for (ParticleVector::const_iterator it=m_particles.begin(); it!=m_particles.end(); it++){
	  const ScanMatcherMap& m1(it->map);
	  const HierarchicalArray2D<PointAccumulator>& h1(m1.storage());
 	  for (int x=0; x<h1.getXSize(); x++){
	    for (int y=0; y<h1.getYSize(); y++){
	      const autoptr< Array2D<PointAccumulator> >& a1(h1.m_cells[x][y]);
	      if (a1.m_reference){
		PointerMap::iterator f=pmap.find(a1.m_reference);
		if (f==pmap.end())
		  pmap.insert(make_pair(a1.m_reference, 1));
		else
		  f->second++;
	      }
	    }
	  }
	}
	cerr << __PRETTY_FUNCTION__ <<  ": Number of allocated chunks" << pmap.size() << endl;
	for(PointerMap::const_iterator it=pmap.begin(); it!=pmap.end(); it++)
	  assert(it->first->shares==(unsigned int)it->second);

	cerr << __PRETTY_FUNCTION__ <<  ": SUCCESS, the error is somewhere else" << endl;
# endif
	GridSlamProcessor* cloned=new GridSlamProcessor(*this);
	
# ifdef MAP_CONSISTENCY_CHECK
	cerr << __PRETTY_FUNCTION__ <<  ": trajectories end" << endl;
	cerr << __PRETTY_FUNCTION__ << ": performing afterclone_fit_test" << endl;
	ParticleVector::const_iterator jt=cloned->m_particles.begin();
	for (ParticleVector::const_iterator it=m_particles.begin(); it!=m_particles.end(); it++){
	  const ScanMatcherMap& m1(it->map);
	  const ScanMatcherMap& m2(jt->map);
	  const HierarchicalArray2D<PointAccumulator>& h1(m1.storage());
	  const HierarchicalArray2D<PointAccumulator>& h2(m2.storage());
	  jt++;
 	  for (int x=0; x<h1.getXSize(); x++){
	    for (int y=0; y<h1.getYSize(); y++){
	      const autoptr< Array2D<PointAccumulator> >& a1(h1.m_cells[x][y]);
	      const autoptr< Array2D<PointAccumulator> >& a2(h2.m_cells[x][y]);
	      assert(a1.m_reference==a2.m_reference);
	      assert((!a1.m_reference) || !(a1.m_reference->shares%2));
	    }
	  }
	}
	cerr << __PRETTY_FUNCTION__ <<  ": SUCCESS, the error is somewhere else" << endl;
# endif
	return cloned;
}
  
  GridSlamProcessor::~GridSlamProcessor(){
    cerr << __PRETTY_FUNCTION__ << ": Start" << endl;
    cerr << __PRETTY_FUNCTION__ << ": Deleting tree" << endl;
    for (std::vector<Particle>::iterator it=m_particles.begin(); it!=m_particles.end(); it++){
#ifdef TREE_CONSISTENCY_CHECK		
      TNode* node=it->node;
      while(node)
	node=node->parent;
      cerr << "@" << endl;
#endif
      if (it->node)
	delete it->node;
      //cout << "l=" << it->weight<< endl;
    }
    
# ifdef MAP_CONSISTENCY_CHECK
    cerr << __PRETTY_FUNCTION__ << ": performing predestruction_fit_test" << endl;
    typedef std::map<autoptr< Array2D<PointAccumulator> >::reference* const, int> PointerMap;
    PointerMap pmap;
    for (ParticleVector::const_iterator it=m_particles.begin(); it!=m_particles.end(); it++){
      const ScanMatcherMap& m1(it->map);
      const HierarchicalArray2D<PointAccumulator>& h1(m1.storage());
      for (int x=0; x<h1.getXSize(); x++){
	for (int y=0; y<h1.getYSize(); y++){
	  const autoptr< Array2D<PointAccumulator> >& a1(h1.m_cells[x][y]);
	  if (a1.m_reference){
	    PointerMap::iterator f=pmap.find(a1.m_reference);
	    if (f==pmap.end())
	      pmap.insert(make_pair(a1.m_reference, 1));
	    else
	      f->second++;
	  }
	}
      }
    }
    cerr << __PRETTY_FUNCTION__ << ": Number of allocated chunks" << pmap.size() << endl;
    for(PointerMap::const_iterator it=pmap.begin(); it!=pmap.end(); it++)
      assert(it->first->shares>=(unsigned int)it->second);
    cerr << __PRETTY_FUNCTION__ << ": SUCCESS, the error is somewhere else" << endl;
# endif
  }


		
  void GridSlamProcessor::setMatchingParameters (double urange, double range, double sigma, int kernsize, double lopt, double aopt, 
						 int iterations, double likelihoodSigma, double likelihoodGain, unsigned int likelihoodSkip){
    m_obsSigmaGain=likelihoodGain;
    m_matcher.setMatchingParameters(urange, range, sigma, kernsize, lopt, aopt, iterations, likelihoodSigma, likelihoodSkip);
    if (m_infoStream)
      m_infoStream << " -maxUrange "<< urange
		   << " -maxUrange "<< range
		   << " -sigma     "<< sigma
		   << " -kernelSize "<< kernsize
		   << " -lstep "    << lopt
		   << " -lobsGain " << m_obsSigmaGain
		   << " -astep "    << aopt << endl;
    
    
  }
  
void GridSlamProcessor::setMotionModelParameters
(double srr, double srt, double str, double stt){
  m_motionModel.srr=srr;
  m_motionModel.srt=srt;
  m_motionModel.str=str;
  m_motionModel.stt=stt;	
  
  if (m_infoStream)
    m_infoStream << " -srr "<< srr 	<< " -srt "<< srt  
		 << " -str "<< str 	<< " -stt "<< stt << endl;
  
}
  
  void GridSlamProcessor::setUpdateDistances(double linear, double angular, double resampleThreshold){
    m_linearThresholdDistance=linear; 
    m_angularThresholdDistance=angular;
    m_resampleThreshold=resampleThreshold;	
    if (m_infoStream)
      m_infoStream << " -linearUpdate " << linear
		   << " -angularUpdate "<< angular
		   << " -resampleThreshold " << m_resampleThreshold << endl;
  }
  
  //HERE STARTS THE BEEF

  GridSlamProcessor::Particle::Particle(const ScanMatcherMap& m):
    map(m), pose(0,0,0), weight(0), weightSum(0), gweight(0), previousIndex(0){
    node=0;
  }
  
  
  void GridSlamProcessor::setSensorMap(const SensorMap& smap){
    
    /*
      Construct the angle table for the sensor
      
      FIXME For now detect the readings of only the front laser, and assume its pose is in the center of the robot 
    */
    
    SensorMap::const_iterator laser_it=smap.find(std::string("FLASER"));
    if (laser_it==smap.end()){
      cerr << "Attempting to load the new carmen log format" << endl;
      laser_it=smap.find(std::string("ROBOTLASER1"));
      assert(laser_it!=smap.end());
    }
    const RangeSensor* rangeSensor=dynamic_cast<const RangeSensor*>((laser_it->second));
    assert(rangeSensor && rangeSensor->beams().size());
    
    m_beams=static_cast<unsigned int>(rangeSensor->beams().size());
    double* angles=new double[rangeSensor->beams().size()];
    for (unsigned int i=0; i<m_beams; i++){
      angles[i]=rangeSensor->beams()[i].pose.theta;
    }
    m_matcher.setLaserParameters(m_beams, angles, rangeSensor->getPose());
    delete [] angles;
  }
  
  void GridSlamProcessor::init(unsigned int size, double xmin, double ymin, double xmax, double ymax, double delta, OrientedPoint initialPose){
    m_xmin=xmin;
    m_ymin=ymin;
    m_xmax=xmax;
    m_ymax=ymax;
    m_delta=delta;
    if (m_infoStream)
      m_infoStream 
	<< " -xmin "<< m_xmin
	<< " -xmax "<< m_xmax
	<< " -ymin "<< m_ymin
	<< " -ymax "<< m_ymax
	<< " -delta "<< m_delta
	<< " -particles "<< size << endl;
    

    m_particles.clear();
    TNode* node=new TNode(initialPose, 0, 0, 0);
    ScanMatcherMap lmap(Point(xmin+xmax, ymin+ymax)*.5, xmax-xmin, ymax-ymin, delta);
    for (unsigned int i=0; i<size; i++){
      m_particles.push_back(Particle(lmap));
      m_particles.back().pose=initialPose;
      m_particles.back().previousPose=initialPose;
      m_particles.back().setWeight(0);
      m_particles.back().previousIndex=0;
      
		// this is not needed
		//		m_particles.back().node=new TNode(initialPose, 0, node, 0);

		// we use the root directly
		m_particles.back().node= node;
    }
    m_neff=(double)size;
    m_count=0;
    m_readingCount=0;
    m_linearDistance=m_angularDistance=0;
  }

  void GridSlamProcessor::processTruePos(const OdometryReading& o){
    const OdometrySensor* os=dynamic_cast<const OdometrySensor*>(o.getSensor());
    if (os && os->isIdeal() && m_outputStream){
      m_outputStream << setiosflags(ios::fixed) << setprecision(3);
      m_outputStream << "SIMULATOR_POS " <<  o.getPose().x << " " << o.getPose().y << " " ;
      m_outputStream << setiosflags(ios::fixed) << setprecision(6) << o.getPose().theta << " " <<  o.getTime() << endl;
    }
  }
   ///////////////////////////:processScan////////////////////////////////////////// 
  /**
   * @brief 
   * 
   * @param reading 
   * @param adaptParticles 
   * @return true 
   * @return false 
   */
  bool GridSlamProcessor::processScan(const RangeReading & reading, int adaptParticles){
    
    /**retireve the position from the reading, and compute the odometry*/
    //得到当前的里程计的位置
    OrientedPoint relPose=reading.getPose();
    //m_count的作用：仅用于区分第0次和其他，第0次是特殊处理,将里程计位姿赋值给激光点云位姿
    if (!m_count){
      m_lastPartPose=m_odoPose=relPose;
    }
    //write the state of the reading and update all the particles using the motion model
    // 对于每一个粒子，都从里程计运动模型中采样,由relPose, m_odoPose的里程计位姿变化量进行更新粒子位姿
    for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++){
      OrientedPoint& pose(it->pose);
    //it->pose  表示机器人上一时刻的最优位置  粒子的位姿(机器人的位姿)
		//relPose   表示里程计算出来的当前位置
		//m_odoPose 表示里程计上一次的位置
    //由relPose, m_odoPose的里程计位姿变化量,加入噪声并进行高斯分布过滤,获得位姿来更新粒子位姿
      pose=m_motionModel.drawFromMotion(it->pose, relPose, m_odoPose);// 运动模型更新------------------------------------------------------------
    }
    //更新文件
    // 根据里程计信息和机器人位姿,update the output file
    if (m_outputStream.is_open()){
      m_outputStream << setiosflags(ios::fixed) << setprecision(6);
      m_outputStream << "ODOM ";
      m_outputStream << setiosflags(ios::fixed) << setprecision(3) << m_odoPose.x << " " << m_odoPose.y << " ";
      m_outputStream << setiosflags(ios::fixed) << setprecision(6) << m_odoPose.theta << " ";
      m_outputStream << reading.getTime();
      m_outputStream << endl;
    }
    if (m_outputStream.is_open()){
      m_outputStream << setiosflags(ios::fixed) << setprecision(6);
      m_outputStream << "ODO_UPDATE "<< m_particles.size() << " ";
      for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++){
	    OrientedPoint& pose(it->pose);
	    m_outputStream << setiosflags(ios::fixed) << setprecision(3) << pose.x << " " << pose.y << " ";
	    m_outputStream << setiosflags(ios::fixed) << setprecision(6) << pose.theta << " " << it-> weight << " ";
      }
      m_outputStream << reading.getTime();
      m_outputStream << endl;
    }
    
    //invoke the callback
    onOdometryUpdate();//没用到
    

    // accumulate the robot translation and rotation
    // 计算机器人平移和旋转增量
    OrientedPoint move=relPose-m_odoPose;
    move.theta=atan2(sin(move.theta), cos(move.theta));
    m_linearDistance+=sqrt(move*move);
    m_angularDistance+=fabs(move.theta);
    
    // if the robot jumps throw a warning
    //如果里程计变化量大于设定阈值,则报警里程计跳变，并发布里程计前后时刻的位姿
    if (m_linearDistance>m_distanceThresholdCheck){
      cerr << "***********************************************************************" << endl;
      cerr << "********** Error: m_distanceThresholdCheck overridden!!!! *************" << endl;
      cerr << "m_distanceThresholdCheck=" << m_distanceThresholdCheck << endl;
      cerr << "Old Odometry Pose= " << m_odoPose.x << " " << m_odoPose.y 
	   << " " <<m_odoPose.theta << endl;
      cerr << "New Odometry Pose (reported from observation)= " << relPose.x << " " << relPose.y 
	   << " " <<relPose.theta << endl;
      cerr << "***********************************************************************" << endl;
      cerr << "** The Odometry has a big jump here. This is probably a bug in the   **" << endl;
      cerr << "** odometry/laser input. We continue now, but the result is probably **" << endl;
      cerr << "** crap or can lead to a core dump since the map doesn't fit.... C&G **" << endl;
      cerr << "***********************************************************************" << endl;
    }
    //更新里程计上一时刻都位姿
    m_odoPose=relPose;
    
    bool processed=false;

    // process a scan only if the robot has traveled a given distance or a certain amount of time has elapsed
    if (! m_count //第一次
       || m_linearDistance>=m_linearThresholdDistance //超过阈值距离
       || m_angularDistance>=m_angularThresholdDistance//超过阈值角度
       || (period_ >= 0.0 && (reading.getTime() - last_update_time_) > period_)){//与上一桢时间相隔超过设定周期
          last_update_time_ = reading.getTime();      //更新上一刻时钟

      if (m_outputStream.is_open()){//文件可读
          m_outputStream << setiosflags(ios::fixed) << setprecision(6);
          // 这是格式化输出流
          // ios::fixed操作符定义在<iotream>中,意思是用小数点形式来显示数据，
          // 有效数字位数默认为为6位，你可以通过setprecision(n)来修改显示有效数字的位数
          m_outputStream << "FRAME " <<  m_readingCount;
          m_outputStream << " " << m_linearDistance;
          m_outputStream << " " << m_angularDistance << endl;
      }
      
      if (m_infoStream)
        	m_infoStream << "update frame " <<  m_readingCount << endl
		     << "update ld=" << m_linearDistance << " ad=" << m_angularDistance << endl;
      
      // cerr与cout的主要区分就是，cout输出的信息可以重定向，而cerr只能输出到标准输出（显示器）上,不需要缓存
      cerr << "Laser Pose= " << reading.getPose().x << " " << reading.getPose().y 
	   << " " << reading.getPose().theta << endl;
      
      
      //this is for converting the reading in a scan-matcher feedable form
      assert(reading.size()==m_beams);//读取激光点云空间
      double * plainReading = new double[m_beams];//创建点云对象,
      for(unsigned int i=0; i<m_beams; i++){
	        plainReading[i]=reading[i];//存储点云数据,用于扫描匹配
      }
      m_infoStream << "m_count " << m_count << endl;//第m_count桢的点云 
      //激光数据复制???????
      RangeReading* reading_copy = 
              new RangeReading(reading.size(),
                               &(reading[0]),
                               static_cast<const RangeSensor*>(reading.getSensor()),
                               reading.getTime());

      if (m_count>0){
        //点云数据
	  scanMatch(plainReading);//--------------------------扫描匹配---------------------------
	  if (m_outputStream.is_open()){//如果文件可以打开
	  m_outputStream << "LASER_READING "<< reading.size() << " ";//显示激光点云数量有效位数2位
	  m_outputStream << setiosflags(ios::fixed) << setprecision(2);
	  for (RangeReading::const_iterator b=reading.begin(); b!=reading.end(); b++){//发布点云
	    m_outputStream << *b << " ";//
	  }
    //读取激光数据并显示
	  OrientedPoint p=reading.getPose();
	  m_outputStream << setiosflags(ios::fixed) << setprecision(6);
	  m_outputStream << p.x << " " << p.y << " " << p.theta << " " << reading.getTime()<< endl;
	  m_outputStream << "SM_UPDATE "<< m_particles.size() << " ";
    // 打印粒子位姿
	  for (ParticleVector::const_iterator it=m_particles.begin(); it!=m_particles.end(); it++){
	    const OrientedPoint& pose=it->pose;
	    m_outputStream << setiosflags(ios::fixed) << setprecision(3) <<  pose.x << " " << pose.y << " ";//坐标精度3位效数字
	    m_outputStream << setiosflags(ios::fixed) << setprecision(6) <<  pose.theta << " " << it-> weight << " ";//角度精度6位有效数字
	  }
	  m_outputStream << endl;
	}
	  onScanmatchUpdate();// 没用到
	
		/**
		 * @brief 由于scanmatch之中对粒子的权重进行了更新，
     * 那么这时候 各个粒子的轨迹上的累计权重都需要重新计算更新树权重。
     * 1、粒子权重归一化,并计算neff值；
     * 2、重置树（将节点的累计权重置0，让节点等于自己的父节点）；
     * 3、重新累加每个粒子的权重
		 */
	  updateTreeWeights(false);
				
      if (m_infoStream){
        m_infoStream << "neff= " << m_neff  << endl;
      }
      if (m_outputStream.is_open()){
        m_outputStream << setiosflags(ios::fixed) << setprecision(6);
        m_outputStream << "NEFF " << m_neff << endl;
      }
      //粒子重采样，根据Neff的大小来进行重采样，不但进行了重采样，也对地图进行了更新
      //重采样的过程里面有第二步的优化，
      resample(plainReading, adaptParticles, reading_copy);
	
  } else {//如果文件不可以打开
              m_infoStream << "Registering First Scan"<< endl;
              for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++){	
              m_matcher.invalidateActiveArea();//设为需要进行计算激活区域
              m_matcher.computeActiveArea(it->map, it->pose, plainReading);//计算激活区域，通过激光雷达的数据计算出来哪个地图栅格应该要被更新了。
              m_matcher.registerScan(it->map, it->pose, plainReading);//计算激光束上所有点的熵,如果不更新地图,默认熵为0,只更新hit
              
              // cyr: not needed anymore, particles refer to the root in the beginning!
              TNode* node=new	TNode(it->pose, 0., it->node,  0);
              //node->reading=0;
              node->reading = reading_copy;
              it->node=node;
              
            }
      }
      //		cerr  << "Tree: normalizing, resetting and propagating weights at the end..." ;
      updateTreeWeights(false);
      //		cerr  << ".done!" <<endl;
      
      delete [] plainReading;
      m_lastPartPose=m_odoPose; //update the past pose for the next iteration
      m_linearDistance=0;
      m_angularDistance=0;
      m_count++;
      processed=true;
      
      //keep ready for the next step
      for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++){
	it->previousPose=it->pose;
      }
      
    }
    if (m_outputStream.is_open())
      m_outputStream << flush;
    m_readingCount++;
    return processed;
  }
  ///////////////////////////:processScan//////////////////////////////////////////
  
  std::ofstream& GridSlamProcessor::outputStream(){
    return m_outputStream;
  }
  
  std::ostream& GridSlamProcessor::infoStream(){
    return m_infoStream;
  }
  
  
  int GridSlamProcessor::getBestParticleIndex() const{
    unsigned int bi=0;
    double bw=-std::numeric_limits<double>::max();
    for (unsigned int i=0; i<m_particles.size(); i++)
      if (bw<m_particles[i].weightSum){
	bw=m_particles[i].weightSum;
	bi=i;
      }
    return (int) bi;
  }

  void GridSlamProcessor::onScanmatchUpdate(){}
  void GridSlamProcessor::onResampleUpdate(){}
  void GridSlamProcessor::onOdometryUpdate(){}

  
};// end namespace




