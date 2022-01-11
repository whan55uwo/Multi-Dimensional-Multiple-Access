function [ALL_Seg_Array, ALL_CDF]=cdf_calculate(ALL,Seg_Num)
UP = max(ALL);
LOW = min(ALL);
ALL_Space = (UP-LOW)/(Seg_Num-1);
ALL_Seg_Array = LOW:ALL_Space:UP;
ALL_PDF = histc(ALL,ALL_Seg_Array);%macro
ALL_CDF = cumsum(ALL_PDF)./length(ALL);
ALL_Seg_Array=ALL_Seg_Array';
