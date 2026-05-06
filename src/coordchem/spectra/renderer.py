import numpy as np
import plotly.graph_objects as go

INTENSITY_SCALE = {
    "very strong":  1.00,
    "strong":       0.75,
    "strong broad": 0.70,
    "medium":       0.50,
    "medium broad": 0.45,
    "weak":         0.25,
}

def gaussian(x, center, height, sigma):
    return height * np.exp(-0.5 * ((x - center) / sigma) ** 2)

def build_spectrum(bands, intensities, sigma=20.0, wn_range=(100, 4000)):
    x = np.linspace(wn_range[0], wn_range[1], 3000)
    y = np.zeros_like(x)
    for band, intensity in zip(bands, intensities):
        y += gaussian(x, band.center, intensity, sigma)
    if y.max() > 0:
        y /= y.max()
    return x, y

def plot_spectrum(bands, intensities, title, sigma):
    x, y = build_spectrum(bands, intensities, sigma=sigma)
    fig  = go.Figure()
    fig.add_trace(go.Scatter(
        x=x, y=y,
        mode="lines",
        fill="tozeroy",
        line=dict(color="#2563eb", width=2),
        fillcolor="rgba(37,99,235,0.08)",
    ))
    for band in bands:
        fig.add_vline(
            x=band.center,
            line=dict(color="rgba(100,116,139,0.3)", width=1, dash="dot")
        )
    fig.update_layout(
        title=title,
        xaxis_title="Wavenumber (cm⁻¹)",
        yaxis_title="Relative Intensity",
        xaxis=dict(autorange="reversed", showgrid=True, gridcolor="#f1f5f9"),
        yaxis=dict(range=[0, 1.08], showgrid=True, gridcolor="#f1f5f9"),
        plot_bgcolor="white",
        height=380,
        margin=dict(l=50, r=20, t=50, b=45),
    )
    return fig