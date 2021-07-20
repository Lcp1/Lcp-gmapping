#include "gmapping/sensor/sensor_range/rangesensor.h"

namespace GMapping{

RangeSensor::RangeSensor(std::string name): Sensor(name){}

RangeSensor::RangeSensor(std::string name, unsigned int beams_num, double res, const OrientedPoint& position, double span, double maxrange):Sensor(name),
	m_pose(position), m_beams(beams_num){
	double angle=-.5*res*beams_num;
	for (unsigned int i=0; i<beams_num; i++, angle+=res){
		RangeSensor::Beam& beam(m_beams[i]);
		beam.span=span;//跨度
		beam.pose.x=0;
		beam.pose.y=0;
		beam.pose.theta=angle;//角度
		beam.maxRange=maxrange;//最大长度
	}
	newFormat=0;
	updateBeamsLookup();
}

void RangeSensor::updateBeamsLookup(){
	//记录激光点三角函数值,后期直接查表
	for (unsigned int i=0; i<m_beams.size(); i++){
		RangeSensor::Beam& beam(m_beams[i]);
		beam.s=sin(m_beams[i].pose.theta);
		beam.c=cos(m_beams[i].pose.theta);
	}
}

};
