#ifndef RANGESENSOR_H
#define RANGESENSOR_H

#include <vector>
#include <gmapping/sensor/sensor_base/sensor.h>
#include <gmapping/utils/point.h>

namespace GMapping{
//已经包括激光点的信息(坐标角度,量程,三角函数表)
// RangeSensor定义了两个成员变量，
// 其中m_pose用于记录传感器的位姿，
// m_beams则记录了扫描光束的信息。
// 同时提供了set-get函数用于访问这些成员变量。
class RangeSensor: public Sensor{
	friend class Configuration;
	friend class CarmenConfiguration;
	friend class CarmenWrapper;
	public:
		struct Beam{
			OrientedPoint pose;	//pose relative to the center of the sensor 位置(x,y)+delta
			double span;	//spam=0 indicates a line-like beam             0表示线条的激光点
			double maxRange;	//maximum range of the sensor             最大量程
			double s,c;		//sinus and cosinus of the beam (optimization);角度三角函数
		};	
		RangeSensor(std::string name);
		RangeSensor(std::string name, unsigned int beams, double res, const OrientedPoint& position=OrientedPoint(0,0,0), double span=0, double maxrange=89.0);
		inline const std::vector<Beam>& beams() const {return m_beams;}
		inline std::vector<Beam>& beams() {return m_beams;}
		inline OrientedPoint getPose() const {return m_pose;}
		void updateBeamsLookup();///记录激光点三角函数值,后期直接查表
		bool newFormat;
	protected:
		OrientedPoint m_pose;
		std::vector<Beam> m_beams;
};

};

#endif
