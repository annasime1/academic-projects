function norm=minmax(signal)
norm=(signal-min(signal))/(max(signal)-min(signal));
end