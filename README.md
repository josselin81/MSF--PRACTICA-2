[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=josselin81/MSF--PRACTICA-2)
# Modelado de Sistemas Fisiológicos. Práctica 2: Sistema respiratorio [GONZALEZC21210380]

## Autor
Gonzalez Garcia Josselin

Ingeniería Biomédica, Departamento de Ingeniería Eléctrica y Electrónica, Tecnológico Nacional de México/IT Tijuana. Blvd. Alberto Limón Padilla s/n, Tijuana, C.P. 22454, B.C., México. Email: l21210380@tectijuana.edu.mx

## Objetivo general
Diseñar un controlador que permita formular un protocolo de tratamiento para que un paciente con enfisema (caso) presente la misma presión alveolar que un individuo sano (control).


## Actividades
1. Calcular analíticamente la función de transferencia del sistema pulmonar.
2. Establecer el modelo de ecuaciones integro-diferenciales.
3. Determinar el error en estado estacionario y la estabilidad del sistema en lazo abierto.
4. Construir el diagrama de bloques para formular un protocolo de tratamiento en un paciente (caso) con base en la respuesta de un individuo sano (control).
5. Diseñar el controlador con Simulink utilizando el bloque PID Controller y la herramienta Tune para sintonizar los valores óptimos para cada una de las ganancias kP, kI y kD.
6. Ilustrar el cambio de flujo del aire y el volumen tidal en respuesta de las siguientes formas de onda de presión sinusoidal en la apertura de la vía aérea [P_ao(t)]:
   a) 15 respiraciones por minuto con una amplitud (A) de 2.5 cmH_2O, es decir, respiración normal.
   b) 30 respiraciones por minuto con una amplitud (A) de 1.5 cmH_2O, es decir, respiración elevada o taquipnea.
7. Determinar la respuesta a la función sinusoidal [u(t)=Asin(wt)] en el intervalo t [0,30] (segundos), en Python, Simulink y Multisim en lazo abierto y en lazo cerrado con el controlador.
8. Elaborar el diagrama biológico del sistema con BioRender.com
9. Discutir los resultados obtenidos en la experimentación in silico y elaborar el reporte de la práctica.

## Docente
Dr. Paul A. Valle

Posgrado en Ciencias de la Ingeniería [PCI] y Departamento de Ingeniería Eléctrica y Electrónica [DIEE], Tecnológico Nacional de México/IT Tijuana. Blvd. Alberto Limón Padilla s/n, Tijuana, C.P. 22454, B.C., México. Email: paul.valle@tectijuana.edu.mx

## Lecturas
[1] Paul. A. Valle, Syllabus para la asignatura de Modelado de Sistemas Fisiológicos, Tecnológico Nacional de México/IT Tijuana, Tijuana, B.C., México, 2025. Permalink: https://www.dropbox.com/scl/fi/4gl55ccrjm9yulvziikxs/Modelado-de-Sistemas-Fisiologicos.pdf

[2] M. C. Khoo, Physiological Control Systems Analysis Simulation, and Estimation, 2nd ed. Piscataway, New Jersey, USA: IEEE Press, 2018.
