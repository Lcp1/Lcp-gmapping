#ifndef SENSORREADING_H
#define SENSORREADING_H

#include "gmapping/sensor/sensor_base/sensor.h"
namespace GMapping{
//基类
class SensorReading{
	public:
		SensorReading(const Sensor* s=0, double time=0);//读取传感器时间,和名字
		virtual ~SensorReading();
		inline double getTime() const {return m_time;}//获取时间
		inline void setTime(double t) {m_time=t;}//设置时间
		inline const Sensor* getSensor() const {return m_sensor;}//设置传感器名字
	protected:
		double m_time;
		const Sensor* m_sensor;

};

}; //end namespace
#endif


