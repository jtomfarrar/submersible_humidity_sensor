function d2byte = rd2byte(words)

%6/1/2000

FirstBit  = 32768; %bin2dec('1000000000000000')
Bits2And3 = 24576; %bin2dec('0110000000000000')
Bits2Not3 = 16384; %bin2dec('0100000000000000')
Bits3Not2 = 8192;  %bin2dec('0010000000000000')
Bits4to16 = 8191;  % bin2dec('0001111111111111')

d2byte=NaN*ones(size(words));

Mult=ones(size(words));

Mult(find(bitand(words(:),FirstBit)==FirstBit)) = -1*...
    Mult(find(bitand(words(:),FirstBit)==FirstBit));

Mult(find(bitand(words,Bits2And3) == Bits3Not2 )) = 0.1 * ...
    Mult(find(bitand(words,Bits2And3) == Bits3Not2 )) ;

Mult(find(bitand(words,Bits2And3) == Bits2Not3 )) = 0.01 * ...
    Mult(find(bitand(words,Bits2And3) == Bits2Not3 )) ;

Mult(find(bitand(words,Bits2And3) == Bits2And3 )) = 0.001 * ...
    Mult(find(bitand(words,Bits2And3) == Bits2And3 ));

d2byte = Mult.*bitand(words,Bits4to16);

% disp([num2str(toc,4) ' sec vectorized']);

end