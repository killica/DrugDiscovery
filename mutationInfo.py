class MutationInfo:
    def __init__(self):
        self.atomSwitchMap = {}
        with open('./mutations/atom_switch.txt', 'r') as atomSwitchFile:
            for line in atomSwitchFile:
                line = line.strip()
                parts = line.split(' ')
                if parts[0] not in self.atomSwitchMap:
                    self.atomSwitchMap[parts[0]] = []
                self.atomSwitchMap[parts[0]].append((parts[1], float(parts[2])))

        self.groupSwitchMap = {}
        with open('./mutations/group_switch.txt', 'r') as groupSwitchFile:
            for line in groupSwitchFile:
                line = line.strip()
                parts = line.split()
                if parts[0] not in self.groupSwitchMap:
                    self.groupSwitchMap[parts[0]] = []
                self.groupSwitchMap[parts[0]].append(parts[1])