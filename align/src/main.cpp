// Author: ChuJiani
// Date: 2023-05-17

#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;
using namespace Eigen;

// 常数定义
const double acc_scale = 1.5258789063e-6;     // FSAS
const double gyr_scale = 1.0850694444e-7;     // CPT
const double freq = 100;                      // 采样频率
const double g = 9.7936174;                   // 当地重力加速度
const double w = 7.292115e-5;                 // 地球自转角速度
const double phi = 30.531651244 * M_PI / 180; // 当地纬度

// IMU 数据结构体
struct ImuData {
    double acc_x;
    double acc_y;
    double acc_z;
    double gyro_x;
    double gyro_y;
    double gyro_z;
};

// 字符串分割
vector<string> split(const string &raw, const char delimiter) {
    vector<string> res;

    size_t start = 0;
    size_t end = raw.find(delimiter);

    while (end != string::npos) {
        res.push_back(raw.substr(start, end - start));
        start = end + 1;
        end = raw.find(delimiter, start);
    }
    res.push_back(raw.substr(start, end - start));

    return res;
}

// 数据读取
vector<ImuData> get_imu_data(const string &file_name, const int time_limit) {
    // 初始化列表
    vector<ImuData> imu_datas;

    // 打开文件
    ifstream file(file_name);
    if (!file.is_open()) {
        cerr << "Failed to open file: " << file_name << endl;
    }

    // 逐行处理数据
    string line;
    while (getline(file, line)) {
        // 限制只读取前 time_limit(单位: s) 数据
        if (imu_datas.size() >= int(time_limit * freq)) {
            break;
        }

        // 数据切割
        vector<string> data = split(split(split(line, ';')[1], '*')[0], ',');

        // 数据转换
        ImuData imu_data;
        imu_data.acc_x = stod(data[5]) * acc_scale * freq;
        imu_data.acc_y = stod(data[4]) * acc_scale * freq * -1;
        imu_data.acc_z = stod(data[3]) * acc_scale * freq;
        imu_data.gyro_x = stod(data[8]) * gyr_scale * freq;
        imu_data.gyro_y = stod(data[7]) * gyr_scale * freq * -1;
        imu_data.gyro_z = stod(data[6]) * gyr_scale * freq;

        // 添加到列表
        imu_datas.push_back(imu_data);
    }

    return imu_datas;
}

// 姿态矩阵转欧拉角
Vector3d dcm_to_euler(const Eigen::Matrix3d &Cnb) {
    double phi = atan2(Cnb(2, 1), Cnb(2, 2));
    double theta = asin(-Cnb(2, 0));
    double psi = atan2(Cnb(1, 0), Cnb(0, 0));
    return Vector3d(phi, theta, psi);
}

// 计算姿态
Vector3d calc_init_attitude(const vector<ImuData> &imu_data, Matrix3d Rn) {
    // 统计 IMU 数据的总值与均值
    double imu_data_sum[6] = {0.0}, imu_data_mean[6] = {0.0};
    for (int i = 0; i < imu_data.size(); i++) {
        imu_data_sum[0] += imu_data[i].acc_x;
        imu_data_sum[1] += imu_data[i].acc_y;
        imu_data_sum[2] += imu_data[i].acc_z;
        imu_data_sum[3] += imu_data[i].gyro_x;
        imu_data_sum[4] += imu_data[i].gyro_y;
        imu_data_sum[5] += imu_data[i].gyro_z;
    }
    for (int i = 0; i < 6; i++) {
        imu_data_mean[i] = imu_data_sum[i] / imu_data.size();
    }

    // 初始化几个向量
    Vector3d gb = {-imu_data_mean[1], -imu_data_mean[0], imu_data_mean[2]};
    gb.normalize();
    Vector3d wb = {imu_data_mean[4], imu_data_mean[3], -imu_data_mean[5]};
    wb.normalize();
    Vector3d vb = gb.cross(wb);
    vb.normalize();
    Matrix3d Rb;
    Rb << gb, wb, vb;
    Rb.transposeInPlace();

    // 姿态计算
    Matrix3d Cnb = Rn.inverse() * Rb;
    return dcm_to_euler(Cnb);
}

int main() {
    // 初始化几个向量
    Vector3d gn = {0, 0, g};
    gn.normalize();
    Vector3d wn = {w * cos(phi), 0, -w * sin(phi)};
    wn.normalize();
    Vector3d vn = gn.cross(wn);
    vn.normalize();
    Matrix3d Rn;
    Rn << gn, wn, vn;
    Rn.transposeInPlace();

    // 读取文件数据
    const string file_name = "data/2.ASC";
    vector<ImuData> imu_data = get_imu_data(file_name, 300);

    // 计算整段时间的平均初始姿态
    Vector3d init_attitude = calc_init_attitude(imu_data, Rn);
    cout << "init attitude(degree): " << init_attitude.transpose() * 180 / M_PI << endl;

    return 0;
}