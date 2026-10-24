from pathlib import Path
import os
import sys
import re
import numpy as np
import pyvista as pv

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"


def natural_key(path):
    return [
        int(x) if x.isdigit() else x.lower()
        for x in re.split(r"(\d+)", path.name)
    ]


files = sorted(
    RESULTS.glob("reentrant*.vtu"),
    key=natural_key
)

if not files:
    raise FileNotFoundError(
        "No VTU files found in results/"
    )

print(f"Found {len(files)} VTU files")

first = pv.read(files[0])

print("Point arrays:", list(first.point_data.keys()))
print("Cell arrays:", list(first.cell_data.keys()))


def find_array(mesh, candidates, components=None):
    for name in candidates:
        for actual in mesh.array_names:
            if actual.lower() == name.lower():
                arr = np.asarray(mesh[actual])
                if components is None or (
                    arr.ndim == 2
                    and arr.shape[1] == components
                ):
                    return actual
    return None


def displacement_name(mesh):
    return find_array(
        mesh,
        ["U", "DISP", "Displacement", "Displacements"],
        components=3
    )


def stress_name(mesh):
    return find_array(
        mesh,
        ["Mises", "von Mises", "S_Mises", "S: Mises",
         "S-Mises", "S.Mises", "S_MISES", "VM"],
    )


def make_frame(index):
    mesh = pv.read(files[index])

    u_name = displacement_name(mesh)
    s_name = stress_name(mesh)

    if u_name:
        deformed = mesh.warp_by_vector(
            u_name,
            factor=1.0
        )
    else:
        deformed = mesh.copy()

    return deformed, u_name, s_name


# Cloud containers typically have no X server, EGL or OSMesa.
# Do not construct any VTK render window in that environment: VTK may
# segfault in native code before Python can catch an exception.
if sys.platform.startswith("linux") and not os.environ.get("DISPLAY"):
    print("Headless environment detected: no DISPLAY.")
    print("Skipping VTK rendering (no X/EGL/OSMesa available).")
    print("VTU files are already ParaView-compatible; open them locally.")
    print("VTU directory:", files[0].parent.resolve())
    for index in sorted({0, len(files)//2, len(files)-1}):
        mesh, u_name, s_name = make_frame(index)
        print(f"Frame {index+1}/{len(files)}: points={mesh.n_points}, "
              f"cells={mesh.n_cells}, displacement={u_name}, stress={s_name}")
    sys.exit(0)

plotter = pv.Plotter()
plotter.set_background("white")

state = {"actor": None}


def update_frame(value):
    index = int(round(value))

    mesh, u_name, s_name = make_frame(index)

    if state["actor"] is not None:
        plotter.remove_actor(state["actor"])

    if s_name:
        state["actor"] = plotter.add_mesh(
            mesh,
            scalars=s_name,
            show_edges=False,
            cmap="turbo",
            scalar_bar_args={
                "title": "Von Mises stress"
            }
        )
    else:
        state["actor"] = plotter.add_mesh(
            mesh,
            color="silver",
            show_edges=True,
            edge_color="gray"
        )

    plotter.add_text(
        f"Frame {index + 1}/{len(files)}\n"
        f"Displacement: {u_name or 'NOT FOUND'}\n"
        f"Stress: {s_name or 'NOT FOUND'}",
        name="frame_info",
        position="upper_left",
        font_size=10,
        color="black"
    )

    plotter.render()


update_frame(0)

if len(files) > 1:
    plotter.add_slider_widget(
        update_frame,
        rng=[0, len(files) - 1],
        value=0,
        title="FEA increment",
        pointa=(0.15, 0.08),
        pointb=(0.85, 0.08),
        style="modern"
    )

plotter.add_axes()
plotter.view_isometric()
plotter.show()
