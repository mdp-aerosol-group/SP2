module SP2

using PyCall
using Dates

function __init__()
    py"""
    import xarray as xr
    import struct
    import numpy as np
    import platform

    from datetime import datetime

    def read_sp2(file_name, debug=False):
        my_data = open(file_name, "rb").read()
        split_file_name = file_name.split("/")
        dt = datetime.strptime(split_file_name[-1][0:8], "%Y%m%d")

        if len(my_data) > 0:
            bytepos = 0
            numCols = struct.unpack(">I", my_data[bytepos:bytepos + 4])[0]
            bytepos += 4
            numChannels = struct.unpack(">I", my_data[bytepos:bytepos + 4])[0]
            if debug:
                print(("Loaded file with numCols = {}, numChannels = {}"
                       .format(numCols, numChannels)))

            data_points_per_record = numChannels * numCols

            bytes_per_record = 2 * data_points_per_record
            bytes_not_data_array = 12 + 2 + 28 + 16
            bytes_per_record += bytes_not_data_array
            last_pos = int(bytes_per_record - 1)
            num_spare_cols = struct.unpack(">I", my_data[last_pos - 4:last_pos])[0]
            if debug:
                print("Number of spare columns = %d" % num_spare_cols)

            if num_spare_cols != 0:
                bytes_per_record += num_spare_cols

            numRecords = int(len(my_data) / bytes_per_record)
            totalRows = numChannels * numRecords
            DataWave = np.zeros((totalRows, numCols), dtype='int16')
            Flag = np.zeros(int(totalRows / numChannels), dtype='int16')
            TimeWave = np.zeros(numRecords, dtype='float64')
            Res1 = np.zeros(numRecords, dtype='float32')
            EventIndex = np.zeros(numRecords, dtype='float32')
            TimeDiv10000 = np.zeros(numRecords, dtype='float64')
            TimeRemainder = np.zeros(numRecords, dtype='float64')
            Res5 = np.zeros(numRecords, dtype='float32')
            Res6 = np.zeros(numRecords, dtype='float32')
            Res7 = np.zeros(numRecords, dtype='float64')
            Res8 = np.zeros(numRecords, dtype='float64')
            if num_spare_cols != 0:
                SpareDataArray = np.zeros(numRecords, num_spare_cols)

            arrayFmt = ">"
            for i in range(data_points_per_record):
                arrayFmt += "h"

            for record in range(numRecords):
                dataStartPoint = record * bytes_per_record + 8
                startRow = record * numChannels
                endRow = startRow + numChannels - 1
                the_row = np.array(struct.unpack(
                    arrayFmt, my_data[dataStartPoint:dataStartPoint + int(data_points_per_record * 2)]))

                DataWave[startRow:endRow + 1, 0:numCols] = the_row.reshape(
                    numCols, numChannels).T
                dataStartPoint += data_points_per_record * 2
                Flag[record] = struct.unpack(">h", my_data[dataStartPoint:dataStartPoint + 2])[0]
                next_floats = struct.unpack(">ffffffff", my_data[dataStartPoint + 2:dataStartPoint + 34])
                TimeWave[record] = next_floats[0]
                Res1[record] = next_floats[1]
                EventIndex[record] = next_floats[2]
                TimeDiv10000[record] = next_floats[3]
                TimeRemainder[record] = next_floats[4]
                Res5[record] = next_floats[5]
                Res6[record] = next_floats[6]
                next_doubles = struct.unpack(">dd", my_data[dataStartPoint + 34:dataStartPoint + 50])
                Res7[record] = next_doubles[0]
                Res8[record] = next_doubles[1]
                dataStartPoint += 50

                if num_spare_cols != 0:
                    startRow = (2 * num_spare_cols) * record
                    dataStartPoint += bytes_not_data_array - 4
                    spareFmt = ">"
                    for i in range(num_spare_cols):
                        spareFmt += "f"

                    SpareDataArray[record] = np.array(
                        struct.unpack(spareFmt, my_data[dataStartPoint:dataStartPoint+4*num_spare_cols]))

            UTCtime = TimeDiv10000 * 10000 + TimeRemainder
            diff_epoch_1904 = (
                datetime(1970, 1, 1) - datetime(1904, 1, 1)).total_seconds()
            UTCdatetime = np.array([
                datetime.fromtimestamp(x - diff_epoch_1904) for x in UTCtime])

            DateTimeWave = (dt - datetime(1904, 1, 1)).total_seconds() + TimeWave

            # Make an xarray dataset for SP2
            Flag = xr.DataArray(Flag, dims={'event_index': EventIndex})
            Res1 = xr.DataArray(Res1, dims={'event_index': EventIndex})
            Res5 = xr.DataArray(Res5, dims={'event_index': EventIndex})
            Res6 = xr.DataArray(Res6, dims={'event_index': EventIndex})
            Res7 = xr.DataArray(Res7, dims={'event_index': EventIndex})
            Res8 = xr.DataArray(Res8, dims={'event_index': EventIndex})
            Time = xr.DataArray(UTCdatetime, dims={'event_index': EventIndex})
            EventInd = xr.DataArray(EventIndex, dims={'event_index': EventIndex})
            DateTimeWaveUTC = xr.DataArray(UTCtime, dims={'event_index': EventIndex})
            DateTimeWave = xr.DataArray(DateTimeWave, dims={'event_index': EventIndex})
            TimeWave = xr.DataArray(TimeWave, dims={'event_index': EventIndex})
            my_ds = xr.Dataset({'time': Time, 'Flag': Flag, 'Res1': Res1, 'Res5': Res5,
                                'Res6': Res6, 'Res7': Res7, 'Res8': Res8, 'EventIndex': EventInd,
                                'DateTimeWaveUTC': DateTimeWaveUTC, 'TimeWave': TimeWave,
                                'DateTimeWave': DateTimeWave})

            for i in range(numChannels):
                temp_array = np.zeros((numRecords, numCols), dtype='int')
                for j in range(numRecords):
                    k = i + j*numChannels
                    temp_array[j] = DataWave[k]
                my_ds['Data_ch' + str(i)] = xr.DataArray(
                    temp_array, dims={'event_index': EventIndex, 'columns': np.arange(0, 100, 1)})
            del my_data
            del DataWave
            return my_ds
        else:
            return None
    """
end

read_sp2(file) = py"read_sp2"(file)

function labview2julia(t)
    basetime = Dates.value(DateTime(1904, 1, 1, 0, 0, 0))
    return t * 1000.0 |> round |> x -> convert(DateTime, Millisecond(x + basetime))
end

function read(file)
    data = read_sp2(file)

    flag = data.Flag.to_numpy()
    eventIndex = data.EventIndex.to_numpy()
    res1 = data.Res1.to_numpy()
    res5 = data.Res5.to_numpy()
    res6 = data.Res6.to_numpy()
    res7 = data.Res7.to_numpy()
    res8 = data.Res8.to_numpy()

    ch0 = data.Data_ch0.to_numpy()
    ch1 = data.Data_ch1.to_numpy()
    ch2 = data.Data_ch2.to_numpy()
    ch3 = data.Data_ch3.to_numpy()
    ch4 = data.Data_ch4.to_numpy()
    ch5 = data.Data_ch5.to_numpy()
    ch6 = data.Data_ch6.to_numpy()
    ch7 = data.Data_ch7.to_numpy()

    s = data.DateTimeWave.to_index()
    t = map(i -> labview2julia(s[i]), 1:length(s))

    return (
        DateTime = t,
        res1 = res1,
        res5 = res5,
        res6 = res6,
        res7 = res7,
        res8 = res8,
        ch0 = ch0,
        ch1 = ch1,
        ch2 = ch2,
        ch3 = ch3,
        ch4 = ch4,
        ch5 = ch5,
        ch6 = ch6,
        ch7 = ch7,
    )
end

end
