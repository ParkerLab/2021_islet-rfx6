class ChunkedFileReader:
    
    def __init__(self, fh, field):
        assert(isinstance(field, int))
        self.fh = fh
        self.field = field
        self.current_value = None
        self.next_set = []
        self.eof = False


    def __next__(self):
        # read until eof or until the field changes
        lines = self.next_set.copy()
        self.next_set = []
        if self.eof:
            raise StopIteration
        
        while True:
            line = self.fh.readline()
            if line == '':
                self.eof = True
                return lines
            line = line.rstrip()
            new_value = line.split()[self.field]
            if self.current_value is None or new_value == self.current_value:
                self.current_value = new_value
                lines.append(line)
            else:
                self.next_set.append(line)
                self.current_value = new_value
                return lines
    
    def __iter__(self):
        return self
