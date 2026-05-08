from cardioseg.gui.main_viewer import launch_viewer
import napari

def main():
    viewer = launch_viewer()
    napari.run()

if __name__ == "__main__":
    main()