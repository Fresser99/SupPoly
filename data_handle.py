import pandas as pd
import os
from datetime import datetime, timedelta
import glob


def generate_timestamp_series(start_date, end_date, interval_minutes):
    """生成时间戳序列"""
    timestamps = []
    current = start_date
    while current <= end_date:
        timestamps.append(current)
        current += timedelta(minutes=interval_minutes)
    return timestamps


def read_txt_file(file_path):
    """读取txt文件并转换为DataFrame"""
    try:
        # 尝试不同的分隔符
        for separator in ['\t', ' ', ',']:
            try:
                df = pd.read_csv(file_path, sep=separator, engine='python')
                # 如果成功读取并且列数正确（3或4列），返回该DataFrame
                if len(df.columns) in [2, 3, 4]:
                    return df
            except:
                continue
        # 如果上述方法都失败，尝试固定宽度读取
        df = pd.read_fwf(file_path)
        if len(df.columns) in [3, 4]:
            return df
        raise ValueError(f"无法正确读取文件{file_path}")
    except Exception as e:
        print(f"读取文件 {file_path} 时出错: {str(e)}")
        return None


def process_source_file(file_path):
    """处理源数据文件，返回时间戳和值的映射"""
    df = read_txt_file(file_path)

    if df is None:
        return {}

    try:
        # 判断文件列数
        if len(df.columns) == 4:  # 4列格式：名称,时间戳,其他,值
            timestamp_col = df.iloc[:, 1]  # 第二列为时间戳
            value_col = df.iloc[:, -2]  # 最后一列为值
        elif len(df.columns) == 3:  # 3列格式：时间戳,其他,值
            timestamp_col = df.iloc[:, 0]  # 第一列为时间戳
            value_col = df.iloc[:, -2]  # 最后一列为值
        elif len(df.columns) == 2:  # 3列格式：时间戳,其他,值
            timestamp_col = df.iloc[:, 0]  # 第一列为时间戳
            value_col = df.iloc[:, -1]  # 最后一列为值
        else:
            print(f"警告: 文件 {file_path} 的列数不正确")
            return {}

        # 确保时间戳列为datetime格式
        timestamps = pd.to_datetime(timestamp_col)

        # 创建时间戳到值的映射
        return dict(zip(timestamps, value_col))
    except Exception as e:
        print(f"处理文件 {file_path} 时出错: {str(e)}")
        return {}


def main(template_file, source_folder, output_file):
    """主函数"""
    try:
        # 读取CSV模板文件获取列名
        template_df = pd.read_csv(template_file)
        columns = template_df.columns.tolist()

        # 生成时间序列
        start_date = datetime(2024, 1, 1)
        end_date = datetime(2024, 9, 1)
        timestamps = generate_timestamp_series(start_date, end_date, 5)

        # 创建结果DataFrame，初始化时间戳列
        result_df = pd.DataFrame(timestamps, columns=[columns[0]])

        # 处理每个需要填充的列
        for col in columns[1:]:
            # 在源文件夹中查找匹配的文件
            matching_files = glob.glob(os.path.join(source_folder, f"*{col}*.txt"))

            if not matching_files:
                print(f"警告: 未找到列 '{col}' 的匹配文件")
                result_df[col] = None
                continue

            print(f"\n处理列 '{col}':")
            # 合并所有匹配文件的数据
            combined_data = {}
            for file in matching_files:
                print(f"正在处理文件: {os.path.basename(file)}")
                file_data = process_source_file(file)
                print(f"- 从文件中读取了 {len(file_data)} 个数据点")
                # 更新combined_data，如果时间戳已存在则使用新值覆盖
                combined_data.update(file_data)

            # 填充数据到结果DataFrame
            values = [combined_data.get(timestamp, None) for timestamp in timestamps]
            result_df[col] = values
            print(f"列 '{col}' 处理完成，共有 {sum(1 for v in values if v is not None)} 个有效数据点")

        # 保存结果
        result_df.to_csv(output_file, index=False)
        print(f"\n处理完成，结果已保存至: {output_file}")
        print(f"总行数: {len(result_df)}")
        # 打印每列的非空值数量
        for col in result_df.columns:
            non_null_count = result_df[col].count()
            print(f"列 '{col}' 的非空值数量: {non_null_count}")

    except Exception as e:
        print(f"程序执行出错: {str(e)}")


if __name__ == "__main__":
    # 设置文件路径
    template_file = "taglist.csv"  # 模板CSV文件路径
    source_folder = 'tafdir'  # 源数据文件夹路径
    output_file = "result.csv"  # 输出文件路径

    main(template_file, source_folder, output_file)