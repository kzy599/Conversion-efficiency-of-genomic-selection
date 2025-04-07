import plotly.graph_objects as go
import pandas as pd
import argparse
def process_data(ITYPE,VTYPE,VSE,plotname,ZTYPE):
    data = pd.read_csv("allzt.csv")
    data.loc[(data['denSNP']==0),'imp']="ped"
    data.loc[(data['denSNP']==1250),'imp']="hd"
    data_filt = data[(data['imp'].isin([ITYPE,"ped","hd"]) ) & (data['nGeneration']==20)]
    x_den = data_filt["denSNP"]
    y_ref = data_filt['nRef']
    z_gg = data_filt[VTYPE]
    z_error = data_filt[VSE]
    x = sorted(data_filt['denSNP'].unique())
    y = sorted(data_filt['nRef'].unique())
    z_matrix = data_filt.pivot(index='nRef', columns='denSNP', values=VTYPE).values
    scatter = go.Scatter3d(
        x=x_den,
        y=y_ref,
        z=z_gg,
        mode='markers',
        marker=dict(size=1, color='red'),
        error_z=dict(
            type='data',
            array=z_error,
            visible=True
        )
    )   

    fig = go.Figure(data=[go.Surface(z=z_matrix, x=x, y=y),scatter])
    fig.update_layout(title='3D Plot of'+' ' + ZTYPE, scene=dict(
                    xaxis_title='Panel Density',
                    yaxis=dict(title='Reference Size',tickvals=y),
                    zaxis_title=ZTYPE),
        margin=dict(l=100, r=100, b=0, t=40),
        width = 800,
        height = 600)
    fig.write_image(plotname)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some strings.')
    parser.add_argument('--ITYPE', type=str, help='Bool of Impuation')
    parser.add_argument('--VTYPE', type=str, help='Colnames of the value')
    parser.add_argument('--VSE', type=str, help='Colnames of the value error')
    parser.add_argument('--plotname', type=str, help='The name of the out plot')
    parser.add_argument('--ZTYPE', type=str, help='The name of the z axis')

    args = parser.parse_args()

    # 调用函数并输出结果到文件
    process_data(args.ITYPE, args.VTYPE, args.VSE, args.plotname,args.ZTYPE)
