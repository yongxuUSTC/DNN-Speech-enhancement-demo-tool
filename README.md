<script src="http://api.html5media.info/1.1.8/html5media.min.js"></script>
This code is only for research. If you want to use it for commercial reason, please contact me: yong.xu.ustc@gmail.com

# DNN-Speech-enhancement-demo-tool
Universal Deep neural network based speech enhancement demo and tools, well pre-trained DNN model

the model/weights can be downloaded here: https://drive.google.com/file/d/1Wg14r0m41kWv2ja-stBkq2uZgjOP5yME/view?usp=sharing  or https://drive.google.com/file/d/0B5r5bvRpQ5DRR1lIV1hpZ0RLQ0E/view?usp=sharing 
also in the config/*   , or you can download from  (@ Baidu Yun) http://pan.baidu.com/s/1eRJGrx4

GPU code for Deep neural network (DNN) based speech enhancement can be found here: https://github.com/yongxuUSTC/DNN-for-speech-enhancement

What kinds of noisy speech can the DNN-enh tool enhance ?

It can enhance any kinds of noisy speech, even the real-world noisy speech. one real-world noisy speech enh demo for the movie : http://staff.ustc.edu.cn/~jundu/The%20team/yongxu/demo/IS15.html

The model is trained only on TIMIT data, so it can get the best performance on the TIMIT test set.

The model can get the best performance on English dataset because TIMIT is US-English. But this tool still can be used to enhance the noisy speech in other languages, like Chinese.

You can use multi-language data to retrain this model to get a general DNN-enh tool.

What else can this code use for?

It is designed for any regression tasks, like speech enhancement, ideal binary/ratio mask (IBM/IRM) estimation, audio/music tagging, acoustic event detection, etc.
Please cite the following papers if you use this code:

[1]A Regression Approach to Speech Enhancement Based on Deep Neural Networks.YongXu,JunDu,Li-Rong Dai and Chin-Hui Lee, IEEE/ACM Transactions on Audio,Speech, and Language Processing,P.7-19,Vol.23,No.1, 2015 (2018 IEEE SPS Best paper award, citations > 350)

[2]An Experimental Study on Speech Enhancement Based on Deep Neural Networks.YongXu, JunDu, Li-Rong Dai and Chin-Hui Lee,IEEE signal processing letters, p. 65-68,vol.21,no. 1,January 2014

[3] Multi-Objective Learning and Mask-Based Post-Processing for Deep Neural Network Based Speech Enhancement, Yong Xu, Jun Du, Zhen Huang, Li-Rong Dai, Chin-Hui Lee, Interspeech2015

Some DNN based speech enhancemen demos: 

http://staff.ustc.edu.cn/~jundu/The%20team/yongxu/demo/SE_DNN_taslp.html http://staff.ustc.edu.cn/~jundu/The%20team/yongxu/demo/IS15.html

[<img src="https://img.youtube.com/vi/Z4quxzXR5yw/maxresdefault.jpg" width="50%">](https://youtu.be/Z4quxzXR5yw) <video src="video.mp4" width="320" height="200" controls preload></video>

<audio src="audio.mp3" controls preload></audio> <audio src="audio.mp3" controls preload></audio>
<audio src="audio.mp3" controls preload></audio> <audio src="audio.wav" controls preload></audio>
<audio src="audio.mp3" controls preload></audio> <audio src="audio.wav" controls preload></audio>
