#ifndef RANGEREADING_H
#define RANGEREADING_H

#include <vector>
#include <gmapping/sensor/sensor_base/sensorreading.h>
#include "gmapping/sensor/sensor_range/rangesensor.h"

namespace GMapping{
/**
 * @brief 
 * SensorReading 读取名字
 * 
 */
class RangeReading: public SensorReading, public std::vector<double>{
	public:
		RangeReading(const RangeSensor* rs, double time=0);
		/**
		 * @brief 
		 * n_beams 激光点数量
		 * d 激光点长度
		 * rs 激光点 坐标角度 量程 三角函数值
		 * time 时间
		 */
		RangeReading(unsigned int n_beams, const double* d, const RangeSensor* rs, double time=0);
		virtual ~RangeReading();
		inline const OrientedPoint& getPose() const {return m_pose;}
		inline void setPose(const OrientedPoint& pose) {m_pose=pose;}
		unsigned int rawView(double* v, double density=0.) const;
		std::vector<Point> cartesianForm(double maxRange=1e6) const;
		unsigned int activeBeams(double density=0.) const;
	protected:
		OrientedPoint m_pose;
};

};

#endif
