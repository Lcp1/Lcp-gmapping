#ifndef SENSOR_H
#define SENSOR_H

#include <string>
#include <map>

namespace GMapping{
//设置传感器名字 基类
class Sensor{
	public:
		Sensor(const std::string& name="");//将名字设为全局变量
		virtual ~Sensor();
		inline std::string getName() const {return m_name;}//获取名字
		inline void setName(const std::string& name) {m_name=name;}//设置名字
	protected:
		std::string m_name;
};

typedef std::map<std::string, Sensor*> SensorMap;

}; //end namespace

#endif

