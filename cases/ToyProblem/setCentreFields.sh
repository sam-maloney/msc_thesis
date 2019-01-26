#!/bin/bash
LOW=0.4
HIGH=0.9
funkySetFields -time 0 -keepPatches -field rho -expression "1" -condition "pos().y>${LOW} && pos().y<${HIGH}"
funkySetFields -time 0 -keepPatches -field rho0 -expression "1" -condition "pos().y>${LOW} && pos().y<${HIGH}"
funkySetFields -time 0 -keepPatches -field p -expression "0.75" -condition "pos().y>${LOW} && pos().y<${HIGH}"
funkySetFields -time 0 -keepPatches -field p0 -expression "0.75" -condition "pos().y>${LOW} && pos().y<${HIGH}"
funkySetFields -time 0 -keepPatches -field T -expression "1.05" -condition "pos().y>${LOW} && pos().y<${HIGH}"
funkySetFields -time 0 -keepPatches -field T0 -expression "1.05" -condition "pos().y>${LOW} && pos().y<${HIGH}"
funkySetFields -time 0 -keepPatches -field U -expression "vector(0,-0.394688351072542,0)" -condition "pos().y>${LOW} && pos().y<${HIGH}"
funkySetFields -time 0 -keepPatches -field U0 -expression "vector(0,-0.394688351072542,0)" -condition "pos().y>${LOW} && pos().y<${HIGH}"
