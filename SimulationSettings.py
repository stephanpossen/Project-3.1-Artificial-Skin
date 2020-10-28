

class SimulationSettings:

    def __init__(self):
        self._sensorList = list     # list of all sensors and their position
        self._sensorFiction = None
        self._sensorYoungM = None
        self._sheetSize = None
        self._sheetPos = None
        self._sheetMass = None
        self._sheetSubdiv = None
        self._sheetYoungM = None
        self._sheetPoisson = None
        # TODO complete the list of all settings


    def read(self, path: str):
        with open(path, 'r') as reader:
            line = reader.readline()
            while line != '':
                if line == 'Sensor':
                    pass
                # TODO add sensors to the sensorlist, with their positions
                if line.startswith('sensorFriction'):
                    data = line.split('=')
                    self._sensorFiction = float(data[1])
                if line.startswith('sensorYoungM'):
                    data = line.split('=')
                    self._sensorYoungM = float(data[1])
                if line.startswith('sheetSize'):
                    data = line.split('=')
                    self._sheetSize = float(data[1])
                if line.startswith('sheetPos'):
                    data = line.split('=')
                    self._sheetPos = float(data[1])
                if line.startswith('sheetMass'):
                    data = line.split('=')
                    self._sheetMass = float(data[1])
                if line.startswith('sheetSubdiv'):
                    data = line.split('=')
                    self._sheetSubdiv = float(data[1])
                if line.startswith('sheetYoungM'):
                    data = line.split('=')
                    self._sheetYoungM = float(data[1])
                if line.startswith('sheetPoisson'):
                    data = line.split('=')
                    self._sheetPoisson = float(data[1])

                line = reader.readline()
