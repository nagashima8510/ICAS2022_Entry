#国内滑走路情報読み込み
import pandas as pd
file_name = "RUNWAY_DB5.csv"

def read_runway():
    rw_table = pd.read_csv(file_name, dtype=str)
    
    #確認
    #print(rw_table)
    
    #型変換
    rw_table['TRUE BRG(deg)'] = rw_table['TRUE BRG(deg)'].astype(float)
    rw_table['Length(m)'] = rw_table['Length(m)'].astype(float)
    rw_table['Width(m)'] = rw_table['Width(m)'].astype(float)
    rw_table['Latitude(deg)'] = rw_table['Latitude(deg)'].astype(float)
    rw_table['Longitude(deg)'] = rw_table['Longitude(deg)'].astype(float)
    rw_table['Altitude'] = rw_table['Altitude'].astype(float)

    #X(km)とY(km)は、大圏経路座標系で計算済み
    rw_table['X(km)'] = rw_table['X(km)'].astype(float)
    rw_table['Y(km)'] = rw_table['Y(km)'].astype(float)

    #滑走路の絞り込み
    rw_table[rw_table['Length(m)'] >= 2500]

    #確認
    #print(rw_table)
    return rw_table