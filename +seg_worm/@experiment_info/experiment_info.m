classdef experiment_info < handle
    %
    %   Class:
    %   seg_worm.experiment_info
    %
    %   I think I had wanted this class to be responsible for creating
    %   the 'info' structure ...
    
    properties
    end
    
    methods
    end
    
end

%{
<configuration>
    <info>
        <program>
            <version>2.0.3.1</version>
        </program>
        <stage>
            <type>Zaber</type>
            <steps>
                <equivalent>
                    <microns>
                        <x>20.9973</x>
                        <y>20.9973</y>
                    </microns>
                    <pixels>
                    <x>-82.89074142959834</x>
<y>-82.82917279011834</y>
</pixels>
</equivalent>
</steps>
</stage>
<camera>
<display>
<id>vfw:Microsoft WDM Image Capture (Win32):0</id>
<resolution>
<width>640</width>
<height>480</height>
</resolution>
<frame>
<rate>30.0</rate>
</frame>
</display>
<recording/>
<effects>
<grayscale>
<on>true</on>
<red>0.3</red>
<green>0.59</green>
<blue>0.11</blue>
</grayscale>
<vignette>
<width>640</width>
<height>480</height>
<raster>
<file>
C:\\Users\\Ithai\\CrossModal\\plasticity genes\\off food fed\\20100421\\mec-4 (u253) off food x_2010_04_21__17_19_20__1.info.xml.vignette.dat
</file>
</raster>
<on>true</on>
</vignette>
</effects>
</camera>
<tracker>
<algorithm>
<rate>1</rate>
<delay>275</delay>
<threshold>
<manual>96</manual>
<auto>
<on>false</on>
<deviation>1.5</deviation>
</auto>
<movement>
<pixel>32</pixel>
<size>0</size>
</movement>
</threshold>
</algorithm>
<boundary>
<type>CENTROID</type>
<centroid>
<x>40.0</x>
<y>40.0</y>
</centroid>
<movement>
<x>1.0</x>
<y>1.0</y>
</movement>
<threshold>
<switch>25.0</switch>
</threshold>
</boundary>
<log/>
</tracker>
</info>
</configuration>
%}