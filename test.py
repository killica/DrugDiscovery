import sys
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QSlider, QLabel
from PyQt5.QtCore import Qt

class MainWindow(QWidget):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Slider Example")
        self.setGeometry(100, 100, 300, 200)

        # Create a vertical layout for the window
        layout = QVBoxLayout()

        # Create the slider, set its range, and orientation (horizontal)
        self.slider = QSlider(Qt.Horizontal)
        self.slider.setRange(0, 100)  # Set the range of the slider (0 to 100)
        self.slider.setTickInterval(1)  # Set the tick interval (1 unit between ticks)

        # Create a label to display the current value of the slider
        self.value_label = QLabel("Value: 0")
        
        # Connect the slider's valueChanged signal to a method
        self.slider.valueChanged.connect(self.update_label)

        # Add the slider and label to the layout
        layout.addWidget(self.slider)
        layout.addWidget(self.value_label)

        # Set the layout of the window
        self.setLayout(layout)

    def update_label(self, value):
        # Update the label with the current value of the slider
        self.value_label.setText(f"Value: {value/100}")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())