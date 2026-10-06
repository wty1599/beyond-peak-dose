import io, zlib, zipfile
STAGED_SIZE = 8161018961
STAGED_SHA256 = 'f2056bac7eb4033e2cc6a5b365f79119cd2ac137b5a54c345065bf27105bd9a5'
MEMBER = 'numericitems.csv'
CSV_BYTES = 80177836504
CSV_CRC32 = 1196767734
CSV_COMPRESSED_BYTES = 8161018679
MAX_PHYSICAL_LINE = 1024 ** 2
EXPECTED_HEADER = ['admissionid', 'itemid', 'item', 'tag', 'value', 'unitid', 'unit', 'comment', 'measuredat', 'registeredat', 'registeredby', 'updatedat', 'updatedby', 'islabresult', 'fluidout']

class BoundedCRCReader(io.RawIOBase):

    def __init__(self, source: zipfile.ZipExtFile):
        self.source = source
        self.total = 0
        self.crc = 0
        self.current_line_bytes = 0

    def readable(self) -> bool:
        return True

    def readinto(self, buffer: bytearray) -> int:
        if self.total == CSV_BYTES:
            if self.source.read(1):
                raise RuntimeError('CSV exceeds declared uncompressed-byte cap')
            return 0
        data = self.source.read(min(len(buffer), CSV_BYTES - self.total))
        n = len(data)
        buffer[:n] = data
        self.total += n
        self.crc = zlib.crc32(data, self.crc)
        for index, segment in enumerate(data.split(b'\n')):
            if index:
                self.current_line_bytes = 0
            self.current_line_bytes += len(segment)
            if self.current_line_bytes > MAX_PHYSICAL_LINE:
                raise RuntimeError('CSV physical line exceeded 1 MiB')
        return n
