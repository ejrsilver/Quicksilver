#include <stdio.h>
#include <string>
#include <chrono>
#include <fstream>
#include <pwr.h>
#include <cassert>

using namespace std;

class Metrics
{
private:
    long instanceTimestamp;
    string fileName;

public:
    Metrics(const string s)
    {
        time_t now = time(0);
        instanceTimestamp = now;
        fileName = "data/" + s + to_string(instanceTimestamp) + ".txt";

        ofstream fp;
        fp.open(fileName);
        if (!fp)
        {
            printf("Error opening %s", fileName);
            exit;
        }

        fp << "TIME_STAMP,PSTATE,CSTATE,CSTATE_LIMIT,SSTATE,CURRENT,VOLTAGE,POWER,POWER_LIMIT_MIN,POWER_LIMIT_MAX,FREQ,FREQ_LIMIT_MIN,FREQ_LIMIT_MAX,ENERGY,TEMP,OS_ID,THROTTLED_TIME,THROTTLED_COUNT," << endl;

        fp.close();
    }

    void getMetrics(PWR_Obj self)
    {
        // open uniq datafile.
        ofstream fp;
        fp.open(fileName, ios::app);
        if (!fp)
        {
            printf("Error opening %s", fileName);
            exit;
        }

        time_t now = time(0);
        long currentTimeStamp = now;

        // get PWR Metrics Values
        uint64_t p_state, c_state, c_state_limit, s_state, os_id, throttled_time, throttled_count;
        p_state = c_state = c_state_limit = s_state = os_id = throttled_time = throttled_count = -1;

        double current, voltage, power, power_min, power_max, frequency, frequency_min, frequency_max, energy, temp;
        current = voltage = power = power_min = power_max = frequency = frequency_min = frequency_max = energy = temp = -1.0;

        PWR_Time ts;
        PWR_ObjAttrGetValue(self, PWR_ATTR_PSTATE, &p_state, &ts);
        PWR_ObjAttrGetValue(self, PWR_ATTR_CSTATE, &c_state, &ts);
        PWR_ObjAttrGetValue(self, PWR_ATTR_CSTATE_LIMIT, &c_state_limit, &ts);
        PWR_ObjAttrGetValue(self, PWR_ATTR_SSTATE, &s_state, &ts);
        PWR_ObjAttrGetValue(self, PWR_ATTR_CURRENT, &current, &ts);
        PWR_ObjAttrGetValue(self, PWR_ATTR_VOLTAGE, &voltage, &ts);
        PWR_ObjAttrGetValue(self, PWR_ATTR_POWER, &power, &ts);
        PWR_ObjAttrGetValue(self, PWR_ATTR_POWER_LIMIT_MIN, &power_min, &ts);
        PWR_ObjAttrGetValue(self, PWR_ATTR_POWER_LIMIT_MAX, &power_max, &ts);
        PWR_ObjAttrGetValue(self, PWR_ATTR_FREQ, &frequency, &ts);
        PWR_ObjAttrGetValue(self, PWR_ATTR_FREQ_LIMIT_MIN, &frequency_min, &ts);
        PWR_ObjAttrGetValue(self, PWR_ATTR_FREQ_LIMIT_MAX, &frequency_max, &ts);
        PWR_ObjAttrGetValue(self, PWR_ATTR_ENERGY, &energy, &ts);
        PWR_ObjAttrGetValue(self, PWR_ATTR_OS_ID, &os_id, &ts);
        PWR_ObjAttrGetValue(self, PWR_ATTR_THROTTLED_TIME, &throttled_time, &ts);
        PWR_ObjAttrGetValue(self, PWR_ATTR_THROTTLED_COUNT, &throttled_count, &ts);

        fp << currentTimeStamp << "," << p_state << "," << c_state << "," << c_state_limit << "," << s_state << "," << current << "," << voltage << "," << power << "," << power_min << "," << power_max << "," << frequency << "," << frequency_min << "," << frequency_max << "," << energy << "," << temp << "," << os_id << "," << throttled_time << "," << throttled_count << endl;
        fp.close();
    }
};
