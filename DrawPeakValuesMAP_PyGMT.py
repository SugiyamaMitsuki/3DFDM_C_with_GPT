import pandas as pd
import numpy as np
import pygmt
import os
import imageio.v2 as imageio  # pip install "imageio[ffmpeg]"　brew install ffmpeg
import gc
import io
import matplotlib.pyplot as plt
import tempfile
from concurrent.futures import ThreadPoolExecutor
from functools import partial

def create_map(centerLon, centerLat, size,frame, lat, lon, v, direction,time,scaleRangemin,scaleRangemax,makeMaxfig=False):

    print(f"direction: {direction}")
    print(f"frame: {frame}")
    print(f"vmax: {v.max()}")
    print(f"vmin: {v.min()}")
    print(f"lat: {lat.min()} ~ {lat.max()}")
    print(f"lon: {lon.min()} ~ {lon.max()}")
    print()

    Map_left = centerLon-size#lon.min()
    Map_right = centerLon+size#lon.max()
    Map_up = centerLat+size#lat.max()
    Map_down = centerLat-size#lat.min()

    projection = "M10c"
    resolution='15s' # 解像度 度d分m秒s "01d", "30m", "20m", "15m", "10m", "06m", "05m", "04m", "03m", "02m", "01m", "30s", "15s", "03s", or "01s"
    scaleKM = 25 #km

    fig = pygmt.Figure()
    # 地図を描画
    fig.basemap(
            region      = [Map_left, Map_right, Map_down, Map_up],
            projection  = projection,
            frame       = ['WSen', 'xaf', 'yaf'], #タイトルを付ける場合'WSen+t"Seismicity Test"'
            )
    
    # 標高水深      
    grid_data = pygmt.datasets.load_earth_relief(
                    resolution=resolution,
                    region = [Map_left, Map_right, Map_down, Map_up],
                    )
    # 勾配計算
    gradient_data = pygmt.grdgradient(
                    grid      = grid_data,
                    azimuth   = [(Map_down+Map_up)/2,(Map_left+Map_right)/2], # 方位角
                    normalize = 'e0.7')
    fig.grdimage(
                    region   = [Map_left, Map_right, Map_down, Map_up],
                    grid     = grid_data, 
                    shading  = gradient_data,
                    cmap     = 'geo',#gray geo batlow
                    transparency = 50,        # コマンド全体に影響する透明度設定
                )
    #縮尺スケールの設置
    map_scale = f"{Map_right-(Map_right-Map_left)*0.2}/{Map_down+(Map_up-Map_down)*0.10}/{(Map_down+Map_up)/2}/{scaleKM}"
    # https://www.pygmt.org/latest/api/generated/pygmt.Figure.coast.html
    fig.coast(
                    shorelines  = 'thinner,black@40',
                    area_thresh = '100',
                    resolution  = 'f',# 'c', 'l', 'i', 'h', 'f' の順に高くなる
                    map_scale   = map_scale,
                    transparency = 80,        # コマンド全体に影響する透明度設定
                    borders     = "1/0.5p",
                    rivers      = "a/1p,skyblue,solid",
                    water       = "lightblue",#'90@10',   # 海域をすこし暗くする "skyblue",#
                )        
    fig.coast(shorelines=True)

    # カラーマップを描画
    pygmt.makecpt(cmap="jet", series=[scaleRangemin,scaleRangemax, 0.001], continuous=True)
    # https://docs.generic-mapping-tools.org/6.3/cookbook/cpts.html
    fig.plot(
        x=lon,
        y=lat,
        fill=v,
        style="c0.02c",
        cmap=True,
        transparency=80,
    )
    
    if makeMaxfig==False:
    # 時間と方向をタイトルに追加
        fig.text(
            x=Map_left+0.1,
            y=Map_up-0.1,
            text=f"{direction} Time: {time:.3f} s",
            font="Helvetica-Bold",
            justify="LM",
            pen="black",
            angle=0,
        )
    else:
        fig.text(
            x=Map_left+0.1,
            y=Map_up-0.1,
            text=f"MAX : {direction}",
            font="Helvetica-Bold",
            justify="LM",
            pen="black",
            angle=0,
        )
    # カラーバーを追加
    # カラーバーを追加
    fig.colorbar(position="JBC+w2.5i/0.15i+h+o0/0.25i", frame=["x+lVelocity", "y+lm/s","+e"])

    # 画像を保存
    # output = f"frame_{direction}_{frame:04d}.png"
    # print(f"output: {output}")
    # fig.savefig(output)
    # fig.show()

    return fig

def process_file(time_step, dt, centerLon, centerLat, size, scaleRangemin, scaleRangemax, direction, output_filename):
    time = dt * time_step
    file_path = f"./output/{output_filename}_{time_step:05d}.bin"
    try:
        with open(file_path, 'rb') as f:
            data_columns = np.fromfile(f, dtype=np.double)
        n = len(data_columns) // 6
        data = pd.DataFrame({
            'lat': data_columns[0:n],
            'lon': data_columns[n:2*n],
            'depth': data_columns[2*n:3*n],
            direction: data_columns[3*n:4*n]
        })
        figMap = create_map(centerLon, centerLat, size, frame=time_step, lat=data['lat'], lon=data['lon'], v=data[direction], direction=direction, time=time, scaleRangemin=scaleRangemin, scaleRangemax=scaleRangemax)
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp_file:
            figMap.savefig(tmp_file.name, dpi=300)
            frame = imageio.imread(tmp_file.name)
            os.remove(tmp_file.name)
        return data, frame
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return None, None
    
def main():
    inputpath = "./input/config.csv"
    df = pd.read_csv(inputpath, comment='#', header=None)
    lat_row = df[df[0] == 'lat0']
    lat = lat_row.iloc[0, 1] if not lat_row.empty else None
    lon_row = df[df[0] == 'lon0']
    lon = lon_row.iloc[0, 1] if not lon_row.empty else None
    print('Latitude: ', lat)
    print('Longitude: ', lon)
    centerLon = lon
    centerLat = lat
    size = 1.5
    dt = 0.01
    interval = 10
    N = 1000
    scaleRangemin=-0.2
    scaleRangemax=0.2
    fps = interval*2.5
    
    output_filename = 'snapshot'
    for z_slice in [0]:
        for direction in ['vx', 'vy', 'vz']:
            frames = []
            maxDataFrame = pd.DataFrame()
            for time_step in range(0, interval * N, interval):
                data, frame = process_file(time_step, dt, centerLon, centerLat, size, scaleRangemin, scaleRangemax, direction, output_filename)
                if data is not None and frame is not None:
                    if maxDataFrame.empty:
                        maxDataFrame = data.copy()
                        maxDataFrame[direction] = maxDataFrame[direction].abs()
                    else:
                        data[direction] = data[direction].abs()
                        maxDataFrame = maxDataFrame.combine(data, np.maximum)
                    frames.append(frame)
            video_filename = f"./output/MapPV_{output_filename}_{direction}_Z={z_slice}_video.mp4"
            imageio.mimsave(video_filename, frames, fps=fps)
            frames.clear()
            gc.collect()
            print(f"Video: {video_filename} is saved.")
            figMap = create_map(centerLon, centerLat, size,frame=0, lat=maxDataFrame['lat'], lon=maxDataFrame['lon'], v=maxDataFrame[direction], direction=direction,time=0,
                                scaleRangemin=0,
                                scaleRangemax=scaleRangemax,
                                makeMaxfig=True)
            figMap.savefig(f"./output/MapPV_{output_filename}_{direction}_Z={z_slice}_max.png", dpi=300)

if __name__ == "__main__":
    main()

