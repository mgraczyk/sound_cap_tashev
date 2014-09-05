function sAngle = NormAngle(sAngle)

while (sAngle > 180)
    sAngle = sAngle - 360;
end
while (sAngle < -180)
    sAngle = sAngle + 360;
end

