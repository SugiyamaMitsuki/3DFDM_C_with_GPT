import pandas as pd
import pygmt
import os

def create_map(centerLon, centerLat, size,frame, lat, lon, v, title):

    def create_basemap(centerLon, centerLat, size, projection, resolution):
        Map_left = centerLon - size
        Map_right = centerLon + size
        Map_up = centerLat + size
        Map_down = centerLat - size

        fig = pygmt.Figure()
        fig.basemap(region=[Map_left, Map_right, Map_down, Map_up], projection=projection, frame=["WSen", "xaf", "yaf"])
        
        return fig, Map_left, Map_right, Map_down, Map_up


    def create_earth_relief(region, resolution):
        grid_data = pygmt.datasets.load_earth_relief(resolution=resolution, region=region)
        gradient_data = pygmt.grdgradient(grid=grid_data, azimuth=[(region[2] + region[3]) / 2, (region[0] + region[1]) / 2], normalize="e0.7")
        
        return grid_data, gradient_data


    def create_coast(region, map_scale):
        fig.coast(shorelines="thinner,black@40", area_thresh="100", resolution="f", map_scale=map_scale, transparency=80, borders="1/0.5p", rivers="a/1p,skyblue,solid", water="lightblue")
        fig.coast(shorelines=True)


    def create_cpt(vmin, vmax):
        # increment = max((v.max() - v.min()) / 100, 1)
        # if v.min() == v.max():
        #     v_max = v.max() + 1
        # else:
        #     v_max = v.max()
        v_min = 0
        v_max = 5000
        increment = max((v_max- v_min) / 100, 1)
        pygmt.makecpt(cmap="lisbon"#"batlow",# https://docs.generic-mapping-tools.org/6.3/cookbook/cpts.html
                    ,series=[v_min, v_max, increment], continuous=True,reverse=False)

    projection = "M10c"
    resolution = "15s"
    
    fig, Map_left, Map_right, Map_down, Map_up = create_basemap(centerLon, centerLat, size, projection, resolution)

    grid_data, gradient_data = create_earth_relief([Map_left, Map_right, Map_down, Map_up], resolution)

    fig.grdimage(region=[Map_left, Map_right, Map_down, Map_up], grid=grid_data, shading=gradient_data, cmap="geo", transparency=50)

    map_scale = f"{Map_right - (Map_right - Map_left) * 0.2}/{Map_down + (Map_up - Map_down) * 0.10}/{(Map_down + Map_up) / 2}/{25}"
    create_coast([Map_left, Map_right, Map_down, Map_up], map_scale)

    create_cpt(v.min(), v.max())

    fig.plot(x=lon, y=lat, fill=v, style="c0.02c", cmap=True, transparency=50)

    fig.text(x=Map_left + 0.1, y=Map_up - 0.1, text=title, font="Helvetica-Bold", justify="LM", pen="black", angle=0)

    fig.colorbar(position="JBC+w2.5i/0.15i+h+o0/0.25i", frame=["x+lSwaveSpeed", "y+lm/s"])


    output = f"./output/Map_physicalProperty_z_val={frame}.png"
    fig.savefig(output)

   
def main():

    file_path = f"./output/Output_physicalProperty.csv"
    data = pd.read_csv(file_path)

    inputpath = "./input/config.csv"
    df = pd.read_csv(inputpath, comment='#', header=None)
    # データフレームからlatとlonの値を取得する。
    lat_row = df[df[0] == 'lat0']
    lat = lat_row.iloc[0, 1] if not lat_row.empty else None

    lon_row = df[df[0] == 'lon0']
    lon = lon_row.iloc[0, 1] if not lon_row.empty else None

    print('Latitude: ', lat)
    print('Longitude: ', lon)
    # centerLon = 139.75
    # centerLat = 35.5
    centerLon = lon
    centerLat = lat
    size = 1.5
    


    for z_val, group in data.groupby("z"):
        print(f"z: {z_val}")
        figMap = create_map(centerLon, centerLat, size,
            frame=z_val, lat=group['lat'], lon=group['lon'], v=group['SwaveSpeed'], title=f"SwaveSpeed at z={z_val}")
        print()

if __name__ == "__main__":
    main()
