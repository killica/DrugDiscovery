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
        # with open('log.txt', 'w') as file:
        #     for (first, others) in self.atomSwitchMap.items():
        #         file.write(first + ": ")
        #         for (second, prob) in others:
        #             file.write(second + " " + str(prob))

        self.groupSwitchMap = {}
        with open('./mutations/group_switch.txt', 'r') as groupSwitchFile:
            for line in groupSwitchFile:
                line = line.strip()
                parts = line.split()
                if parts[0] not in self.groupSwitchMap:
                    self.groupSwitchMap[parts[0]] = []
                self.groupSwitchMap[parts[0]].append(parts[1])
        # with open('log1.txt', 'w') as file:
        #     for (first, others) in self.groupSwitchMap.items():
        #         file.write(first + ": ")
        #         for second in others:
        #             file.write(second)