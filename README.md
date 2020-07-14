# Alert-signal-detector-and-separator
# Design and Integration of alert signal detector and separator for hearing Aid applications


## Overview
This GitHub repository provides for Alert signal detection and separation on iOS smartphone platforms. The example app provided here is for hearing improvement studies. 

> **Abstract:** Alert signals like sirens and alarms are important as they warn people of precarious situations. This work presents the detection and separation of these acoustically important alert signals to assist the hearing impaired listeners. The proposed method is based on convolutional neural network (CNN) and convolutional-recurrent neural network (CRNN). The developed method consists of two blocks, the detector block, and the separator block. The entire setup is integrated with speech enhancement (SE) algorithms used in a hearing aid device (HAD) signal processing pipeline. The detector recognizes the presence of the alert signal in various noisy environments. The separator block separates the alert signal from the mixture of noisy signals before passing it through SE to ensure minimal or no attenuation of the alert signal. It is implemented on a smartphone as an application that seamlessly works with HADs in real-time. This smartphone assistive setup allows the hearing aid users to know the presence of the alert sounds even when these are out of sight. The algorithm is computationally efficient with a low processing delay. The key contribution of this paper includes the development and integration of alert signal separator block with SE and the realization of the entire setup on a smartphone in real-time. The proposed method is compared with several state-of-the-art techniques through objective measures in various noisy conditions. The experimental analysis demonstrates the effectiveness and practical usefulness of the developed setup in real-world noisy scenarios.

You can find the paper for this GitHub repository: https://ieeexplore.ieee.org/abstract/document/9106356

## Audio-Video Demo

- https://youtu.be/lwZmtfAhZeQ
- https://utdallas.edu/ssprl/hearing-aid-project/speech-enhancement/alert-signal-integration/

## Users Guides

[iOS](https://github.com/ssprl/Alert-signal-detector-and-separator/blob/master/User%20Guide-%20iOS%20Alert%20Signal%20Detection.docx.pdf)

## Requirements 
- iPhone X running iOS 13.0.1
- Tensorflow 
- TensorflowLite
- Firebase

## License and Citation
The codes are licensed under MIT license.

For any utilization of the code content of this repository, one of the following books needs to get cited by the user:

- G. S. Bhat, N. Shankar and I. M. S. Panahi, "Design and Integration of Alert Signal Detector and Separator for Hearing Aid Applications," in IEEE Access, vol. 8, pp. 106296-106309, 2020, doi: 10.1109/ACCESS.2020.2999546.

## Disclaimer
This work was supported in part by the National Institute of the Deafness and Other Communication Disorders (NIDCD) of the National

Institutes of Health (NIH) under Award 1R01DC015430-04. The content is solely the responsibility of the authors and does not necessarily represent the official views of the NIH
