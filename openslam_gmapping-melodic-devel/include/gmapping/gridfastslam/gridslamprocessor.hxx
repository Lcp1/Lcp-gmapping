
#ifdef MACOSX
// This is to overcome a possible bug in Apple's GCC.
#define isnan(x) (x==FP_NAN)
#endif
////////////////////////////////scanMatch//////////////////////////////////////////////
/**Just scan match every single particle.
If the scan matching fails, the particle gets a default likelihood.*/
// plainReading 点云数据                    
inline void GridSlamProcessor::scanMatch(const double* plainReading){
  // sample a new pose from each scan in the reference
  
  double sumScore=0;//总得分初始化
  //对每个粒子都进行优化
  for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++){
    OrientedPoint corrected;
    double score, l, s;
    // corrected 为返回最高的分位姿 it->pose 为机器人位姿 plainReading点云数据
    /**
     * @brief 通过计算机器人附近最高分的位姿,返回位姿
     * corrected 修正后的最高分位姿
     * it->map 地图
     * it->pose 机器人坐标
     * plainReading 点云信息
     */
    score=m_matcher.optimize(corrected, it->map, it->pose, plainReading);//函数主要入口 优化
    //    it->pose=corrected;
    if (score>m_minimumScore){//如果分值大于有效最小阈值,就进行更新位姿
      it->pose=corrected;
    } else {
        if (m_infoStream){
          m_infoStream << "Scan Matching Failed, using odometry. Likelihood=" << l <<std::endl;
          m_infoStream << "lp:" << m_lastPartPose.x << " "  << m_lastPartPose.y << " "<< m_lastPartPose.theta <<std::endl;
          m_infoStream << "op:" << m_odoPose.x << " " << m_odoPose.y << " "<< m_odoPose.theta <<std::endl;
        }
    }
    /**
    *
    *(s, l得分, it->map地图, it->pose机器人位姿, plainReading激光点云数据)

    */ 
   /**
    * 粒子的最优位姿计算了之后，重新计算粒子的权重和似然，optimize得到了最优位姿，likelihoodAndScore算最优位姿的似然 l即是权重
    * s: 分布概率分数
    * l: 似然
    * it->map : 地图
    * it->pose: 机器人位姿
    * plainReading : 激光点云数据
    * 有个疑问, 如果ScanMatcher::score()函数也顺便把l似然计算,是不是就不需要ScanMatcher::likelihoodAndScore()函数了??????????/???????
    * **/
    m_matcher.likelihoodAndScore(s, l, it->map, it->pose, plainReading);
    sumScore+=score;
    it->weight+=l;
    it->weightSum+=l;

    //set up the selective copy of the active area
    //by detaching the areas that will be updated
    m_matcher.invalidateActiveArea();//设定地图动超出领域进行扩展地图标记
    // 计算有效区域，通过激光雷达的数据计算出来哪个地图栅格应该要被更新了。
    // (这里只是计算出来栅格的位置，然后插入地图中,并不对数据进行更新)
    m_matcher.computeActiveArea(it->map, it->pose, plainReading);
  }
  if (m_infoStream)
    m_infoStream << "Average Scan Matching Score=" << sumScore/m_particles.size() << std::endl;	//打印平均分
}
///////////////////////////////////scanMatch///////////////////////////////////////////
inline void GridSlamProcessor::normalize(){
  //normalize the log m_weights
  double gain=1./(m_obsSigmaGain*m_particles.size());
  double lmax= -std::numeric_limits<double>::max();
  //获取最大权重值
  for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++){
    lmax=it->weight>lmax?it->weight:lmax;
  }
  //cout << "!!!!!!!!!!! maxwaight= "<< lmax << endl;
  m_weights.clear();
  double wcum=0;
  m_neff=0;
  // 权重累加
  for (std::vector<Particle>::iterator it = m_particles.begin(); it != m_particles.end(); it++){
    m_weights.push_back( exp(gain*(it->weight-lmax)));
    wcum+=m_weights.back();
    //cout << "l=" << it->weight<< endl;
  }
  
  m_neff=0;
  //每个点权重归一化
  for (std::vector<double>::iterator it=m_weights.begin(); it!=m_weights.end(); it++){
    *it=*it/wcum;
    double w=*it;
    m_neff+=w*w;
  }
  //求出neff
  m_neff=1./m_neff;
  
}
/**
 * @brief 根据指标Neff进行重采样，更新重采样后各个粒子的位姿和地图。
 * 
 * @param plainReading 激光传感器的原始读数，用于更新粒子的地图
 * @param adaptSize 一个配置参数，用于指定重采样的粒子数量 
 * @param reading  打了时间戳的传感器读数，其中还记录了采集数据时传感器的位姿
 * @return true 
 * @return false 
 */
inline bool GridSlamProcessor::resample(const double* plainReading, int adaptSize, const RangeReading* reading){
  
  bool hasResampled = false;//重采样标记位
  // vector {const OrientedPoint& pose, double weight, TNode* parent=0, unsigned int childs=0}
  TNodeVector oldGeneration;//存储粒子状态容器是机器人运动轨迹树的节点类型，
  // vector<GridSlamProcessor::TNode*>  
  // 从轨迹树的根节点开始到每个叶子节点所描述的路径，就是一个粒子所记录的运动轨迹。
  for (unsigned int i=0; i<m_particles.size(); i++){//存储粒子
    oldGeneration.push_back(m_particles[i].node);
  }
  // 衡量粒子权重的相似性 当小于阈值时，执行重采样,neff值越大相似性越高,就不用重采样---------------------------if
  if (m_neff<m_resampleThreshold*m_particles.size()){		
    if (m_infoStream)
      m_infoStream  << "*************RESAMPLE***************" << std::endl;
    
    uniform_resampler<double, double> resampler;//创建重采样器

    //返回一个std::vector<unsigned int>类型的数组，记录了采样后的样本在原始粒子集合的索引。
    m_indexes=resampler.resampleIndexes(m_weights, adaptSize);//------------------已经重采样完毕
    
    if (m_outputStream.is_open()){
      m_outputStream << "RESAMPLE "<< m_indexes.size() << " ";
      for (std::vector<unsigned int>::const_iterator it=m_indexes.begin(); it!=m_indexes.end(); it++){
      	m_outputStream << *it <<  " ";
      }
      m_outputStream << std::endl;
    }
    
    onResampleUpdate();
    //BEGIN: BUILDING TREE
    ParticleVector temp;//临时粒子集合,用于记录采样后的粒子集合,最终它将被整合到建图引擎的粒子集合m_particles中。
    unsigned int j=0;
    // 这是为了删除已重新取样的粒子。
    std::vector<unsigned int> deletedParticles;  		//this is for deleteing the particles which have been resampled away.
    
    //		cerr << "Existing Nodes:" ;
    // 
    for (unsigned int i=0; i<m_indexes.size(); i++){
      //			cerr << " " << m_indexes[i];
      while(j<m_indexes[i]){//j是比m_indexes小的,当j=m_indexes,序号相同就退出循环
            deletedParticles.push_back(j);//记录m_indexes没有的粒子
            j++;
			}
      if (j==m_indexes[i])//相同后,j指向下一个粒子序号
        	j++;
      // ，为各个粒子更新最新位姿和传感器数据，并记录在轨迹树中。
      Particle & p=m_particles[m_indexes[i]];//粒子
      TNode* node=0;
      TNode* oldNode=oldGeneration[m_indexes[i]];//获取旧的粒子
      //			cerr << i << "->" << m_indexes[i] << "B("<<oldNode->childs <<") ";
      node=new	TNode(p.pose, 0, oldNode, 0);
      //node->reading=0;
      node->reading=reading;
      //			cerr << "A("<<node->parent->childs <<") " <<endl;
      
      temp.push_back(p);//存储位置
      temp.back().node=node;//存储位置对应的节点
      temp.back().previousIndex=m_indexes[i];//previousIndex存储指向前一个重采样后的粒子索引（内部使用）
    }
    // 删除for循环遗漏的粒子
    while(j<m_indexes.size()){
      deletedParticles.push_back(j);
      j++;
    }
    //		cerr << endl;
    std::cerr <<  "Deleting Nodes:";
    //删除被重采样虑掉的粒子,释放内存
    for (unsigned int i=0; i<deletedParticles.size(); i++){
      std::cerr <<" " << deletedParticles[i];
      delete m_particles[deletedParticles[i]].node;
      m_particles[deletedParticles[i]].node=0;
    }
    std::cerr  << " Done" <<std::endl;
    
    //END: BUILDING TREE
    std::cerr << "Deleting old particles..." ;
    m_particles.clear();//清理粒子容器
    std::cerr << "Done" << std::endl;
    std::cerr << "Copying Particles and  Registering  scans...";
    for (ParticleVector::iterator it=temp.begin(); it!=temp.end(); it++){
      it->setWeight(0);
      m_matcher.invalidateActiveArea();//设置处理地图标记位
      m_matcher.registerScan(it->map, it->pose, plainReading);//计算激光束上的熵总和
      m_particles.push_back(*it);//保存粒子
    }
    std::cerr  << " Done" <<std::endl;
    hasResampled = true;//标记完成采样
  } else {//不进行重采样,也需要更新各个粒子的轨迹和地图。------------------------------------------------------else
    int index=0;
    std::cerr << "Registering Scans:";
    TNodeVector::iterator node_it=oldGeneration.begin();
    for (ParticleVector::iterator it=m_particles.begin(); it!=m_particles.end(); it++){
      //create a new node in the particle tree and add it to the old tree
      //BEGIN: BUILDING TREE  
      TNode* node=0;
      /**
       * @brief 
       * p it->pose 机器人位姿
       * w  0.0 权重
       * n node_it 父节点
       * c 0 子节点数量
       */
      // TNode(const OrientedPoint& p, double w, TNode* n, unsigned int c)
      node=new TNode(it->pose, 0.0, *node_it, 0);
      
      //TNode 初始值 node->reading=0;
      node->reading=reading;
      it->node=node;//赋值节点

      //END: BUILDING TREE
      m_matcher.invalidateActiveArea();//激活计算活动区域
      m_matcher.registerScan(it->map, it->pose, plainReading);//记录激光束同时,记录熵
      it->previousIndex=index;
      index++;
      node_it++;//父节点递增
      
    }
    std::cerr  << "Done" <<std::endl;
    
  }
  //END: BUILDING TREE
  
  return hasResampled;
}
