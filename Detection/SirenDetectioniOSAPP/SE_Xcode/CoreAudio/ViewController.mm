//
//  ViewController.m
//  CoreAudio
//
//  Created by Shankar, Nikhil on 4/4/17.
//  Copyright © 2017 default. All rights reserved.
//

#import "ViewController.h"
#import "FIR.h"
#import "overlapAdd.hpp"
#include <chrono>

#import "Firebase/Firebase.h"
#import "FirebaseMLCommon/FirebaseMLCommon.h"
#import "FirebaseMLModelInterpreter/FirebaseMLModelInterpreter.h"

using namespace std;
using namespace std::chrono;


#define kOutputBus 0
#define kInputBus 1
#define SHORT2FLOAT 1/32768.0
#define FLOAT2SHORT 32768.0;
#define FRAMESIZE 256
#define SAMPLINGFREQUENCY 16000
#define MAXIT 1000
#define EULER 0.5772156649
#define FPMIN 1.0e-30
#define EPS  1.0e-7f
#define pi 3.1415926535897932384626433832795
#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif
int nFFT=512;
int BUFSIZE = 256;
int on;
int sir_on;
int frameCounter=0;
static float *output_final = (float*)calloc(FRAMESIZE, sizeof(float));
static float *input = (float*)calloc(FRAMESIZE, sizeof(float));
static float *float_output = (float*)calloc(FRAMESIZE, sizeof(float));
//SIren_detection parameters
static float *siren_fea = (float*)calloc(nFFT + 1,sizeof(float));
NSArray *modelOut;
static float *data = (float*)calloc(nFFT + 1, sizeof(float));

//static short *buffer = (short*)calloc(FRAMESIZE/sizeof(short), sizeof(short));
//static short *output = (short*)calloc(FRAMESIZE/sizeof(short), sizeof(short));

overlapAdd OA;

AURenderCallbackStruct callbackStruct;
AudioUnit au;
AudioBuffer tempBuffer;

@interface ViewController ()

@end

@implementation ViewController

@synthesize EnhancedSwitch;
@synthesize volumeslider;
@synthesize SirenSwitch;
//@synthesize NormaLabel;
@synthesize sirenLabel;



static OSStatus playbackCallback(void *inRefCon,AudioUnitRenderActionFlags *ioActionFlags,const AudioTimeStamp *inTimeStamp,UInt32 inBusNumber, UInt32 inNumberFrames, AudioBufferList *ioData)
{

    for (int i=0; i < ioData->mNumberBuffers; i++) {
        AudioBuffer buffer = ioData->mBuffers[i];
        UInt32 size = min(buffer.mDataByteSize, tempBuffer.mDataByteSize);
        memcpy(buffer.mData, tempBuffer.mData, size);
        buffer.mDataByteSize = size;
    }
    return noErr;
}

static OSStatus recordingCallback(void *inRefCon,AudioUnitRenderActionFlags *ioActionFlags,const AudioTimeStamp *inTimeStamp, UInt32 inBusNumber, UInt32 inNumberFrames, AudioBufferList *ioData)
{

    AudioBuffer buffer;
    ViewController* view = (__bridge ViewController *)(inRefCon);

    buffer.mNumberChannels = 1;
    buffer.mDataByteSize = inNumberFrames * 2;
    buffer.mData = malloc( inNumberFrames * 2 );

    // Put buffer in a AudioBufferList
    AudioBufferList bufferList;
    bufferList.mNumberBuffers = 1;
    bufferList.mBuffers[0] = buffer;

    AudioUnitRender(au, ioActionFlags, inTimeStamp,inBusNumber,inNumberFrames,&bufferList);
    [view processAudio:&bufferList];
    // printf("%f\n",buffer);

    return noErr;
}


void RUN_ON_UI_THREAD(dispatch_block_t block)
{
    if ([NSThread isMainThread])
        block();
    else
        dispatch_sync(dispatch_get_main_queue(), block);
}

- (void)didReceiveMemoryWarning {
    [super didReceiveMemoryWarning];
    // Dispose of any resources that can be recreated.
}

-(void) modelRun:(NSMutableData*)input{

    //Initialise Model
    NSString *modelPath = [NSBundle.mainBundle pathForResource:@"CNN_Weights-048"
                                                        ofType:@"tflite"];
    FIRCustomLocalModel *localModel =
    [[FIRCustomLocalModel alloc] initWithModelPath:modelPath];

    FIRModelInterpreter *interpreter =
    [FIRModelInterpreter modelInterpreterForLocalModel:localModel];

    FIRModelInputOutputOptions *ioOptions = [[FIRModelInputOutputOptions alloc] init];
    NSError *error;
    [ioOptions setInputFormatForIndex:0
                                 type:FIRModelElementTypeFloat32
                           dimensions:@[@1, @257]
                                error:&error];
    if (error != nil) { return; }
    [ioOptions setOutputFormatForIndex:0
                                  type:FIRModelElementTypeFloat32
                            dimensions:@[@1, @2]
                                 error:&error];
    if (error != nil) { return; }

    // Initialise Inderence
    FIRModelInputs *inputs = [[FIRModelInputs alloc] init];
    [inputs addInput:input error:&error];

    [interpreter runWithInputs:inputs options:ioOptions completion:^(FIRModelOutputs * _Nullable outputs, NSError * _Nullable error)
     {
//                        if (error != nil || outputs == nil) {
//                            return;
//                        }
                        NSError *outputError;
                        NSArray *probabilites = [outputs outputAtIndex:0 error:&outputError][0];
                        modelOut = probabilites;
                        //NSLog (@ "%@",modelOut[0]);

    }];
}

-(void) updatemodel{
        ViewController *sampleClass = [[ViewController alloc]init];


        NSMutableData *inputData = [[NSMutableData alloc] initWithCapacity:0];
        for (int row = 0; row < 1; row++) {
            for (int col = 0; col < 257; col++) {
                  Float32 red = data[col];
                  //printf( "%f\n",red);
                  [inputData appendBytes:&red length:sizeof(red)];
            }
        }
        [sampleClass modelRun:inputData];
    //NSLog(@ "%@",modelOut[0]);
}


-(void) processAudio: (AudioBufferList*) bufferList{
    AudioBuffer sourceBuffer = bufferList->mBuffers[0];
   // short *buffer = (short*)calloc(sourceBuffer.mDataByteSize);
    //short *output =(short*)calloc(sourceBuffer.mDataByteSize);
   // static float *input1 = (float*)calloc(nFFT, sizeof(float));
   // static float *store_buffer = (float*)calloc(FRAMESIZE, sizeof(float));
    static short *buffer = (short*)calloc(sourceBuffer.mDataByteSize/sizeof(short), sizeof(short));
    static short *output = (short*)calloc(sourceBuffer.mDataByteSize/sizeof(short), sizeof(short));
    //static float *float_output1 = (float*)calloc(nFFT, sizeof(float));

    memcpy(buffer, bufferList->mBuffers[0].mData, bufferList->mBuffers[0].mDataByteSize);

    for (int i = 0; i < FRAMESIZE; i++)
    {
        input[i] = buffer[i] * SHORT2FLOAT;
    }
    if(on==0)
    {
        for(int i =0;i<FRAMESIZE;i++)
        {
            float_output[i]=input[i];
        }
    }
    if(on==1)
    {
        if(frameCounter==0)
        {
            OA.setFrameSize(BUFSIZE);
            OA.setFftSize(nFFT);
            OA.init();
            //[self updatemodel];
        }
        OA.process(input, output_final, frameCounter);
        OA.siren_det(input, siren_fea, frameCounter);
        for (int i = 0; i < nFFT + 1; i++)
        {
            data[i] = siren_fea[i];
        }
//        auto start = high_resolution_clock::now();
        if (sir_on==1){
        [self updatemodel];
        
//        auto stop = high_resolution_clock::now();
        //auto duration = duration_cast<microseconds>(stop - start);
        //printf("Time taken by function: %lld microseconds\n", duration.count() );
        //NSLog(@ "%@",modelOut[0]);
        if(modelOut[0]>modelOut[1])
        {
            RUN_ON_UI_THREAD(^{
                //noiseLabel.text = @"Noise";
                //NormaLabel.backgroundColor = [UIColor greenColor];
                sirenLabel.backgroundColor = [UIColor blackColor];

            });
        }
         else if (modelOut[1]>modelOut[0])
        {
            RUN_ON_UI_THREAD(^{
                //noiseLabel.text = @"Speech";
                sirenLabel.backgroundColor = [UIColor redColor];
                //NormaLabel.backgroundColor = [UIColor blackColor];
            });
        }
        }
        frameCounter++;
        //fir(output_final,float_output,FRAMESIZE);

        for (int i = 0; i < FRAMESIZE; i++)
        {
            float_output[i] = output_final[i];
        }
        ///fir(output_final,float_output,FRAMESIZE);
    }
    for (int i = 0; i < FRAMESIZE; i++)
    {
        float_output[i]=float_output[i]*2;
        output[i] = float_output[i] * FLOAT2SHORT ;
    }
    if (tempBuffer.mDataByteSize != sourceBuffer.mDataByteSize)
    {
        free(tempBuffer.mData);
        tempBuffer.mDataByteSize = sourceBuffer.mDataByteSize;
        tempBuffer.mData = malloc(sourceBuffer.mDataByteSize);
    }
    memcpy(tempBuffer.mData, output, bufferList->mBuffers[0].mDataByteSize);
}

//-(void) updatesnr{
//    snrdisplay.text=[ NSString stringWithFormat:@"%f",y];
//
//}


- (void)viewDidLoad {
    [super viewDidLoad];
    volumeslider.value=[[MPMusicPlayerController applicationMusicPlayer] volume];
    //betaLabel.text = @"0.5";
    [[AVAudioSession sharedInstance] setCategory: AVAudioSessionCategoryPlayAndRecord error: NULL];
    [[AVAudioSession sharedInstance] setMode: AVAudioSessionModeVideoRecording error:NULL];
    [[AVAudioSession sharedInstance] setPreferredSampleRate:SAMPLINGFREQUENCY error:NULL];
    [[AVAudioSession sharedInstance]
     setPreferredIOBufferDuration:(float)FRAMESIZE/(float)SAMPLINGFREQUENCY error:NULL];
    AudioComponentDescription desc;
    desc.componentType = kAudioUnitType_Output;
    desc.componentSubType = kAudioUnitSubType_RemoteIO;
    desc.componentFlags = 0;
    desc.componentFlagsMask = 0;
    desc.componentManufacturer = kAudioUnitManufacturer_Apple;
    AudioComponent component = AudioComponentFindNext(NULL, &desc);
    if (AudioComponentInstanceNew(component, &au) != 0) abort();

    UInt32 value = 1;
    if (AudioUnitSetProperty(au, kAudioOutputUnitProperty_EnableIO, kAudioUnitScope_Output, 0, &value,
                             sizeof(value))) abort();
    value = 1;
    if (AudioUnitSetProperty(au, kAudioOutputUnitProperty_EnableIO, kAudioUnitScope_Input, 1, &value,
                             sizeof(value))) abort();



    AudioStreamBasicDescription format;
    format.mSampleRate = 0;
    format.mFormatID = kAudioFormatLinearPCM;
    format.mFormatFlags = kAudioFormatFlagIsSignedInteger;
    format.mFramesPerPacket = 1;
    format.mChannelsPerFrame = 1;
    format.mBitsPerChannel = 16;
    format.mBytesPerPacket = 2;
    format.mBytesPerFrame = 2;

    if (AudioUnitSetProperty(au, kAudioUnitProperty_StreamFormat, kAudioUnitScope_Input, 0, &format,
                             sizeof(format))) abort();
    if (AudioUnitSetProperty(au, kAudioUnitProperty_StreamFormat, kAudioUnitScope_Output, 1, &format,
                             sizeof(format))) abort();
    // Set input callback

    callbackStruct.inputProc = recordingCallback;
    callbackStruct.inputProcRefCon = (__bridge void *)(self);
    AudioUnitSetProperty(au, kAudioOutputUnitProperty_SetInputCallback, kAudioUnitScope_Global, kInputBus,  &callbackStruct, sizeof(callbackStruct));

    // Set output callback
    callbackStruct.inputProc = playbackCallback;
    callbackStruct.inputProcRefCon = (__bridge void *)(self);
    AudioUnitSetProperty(au, kAudioUnitProperty_SetRenderCallback, kAudioUnitScope_Global, kOutputBus,&callbackStruct, sizeof(callbackStruct));
    tempBuffer.mNumberChannels = 1;
    tempBuffer.mDataByteSize = FRAMESIZE * 2;
    tempBuffer.mData = malloc( FRAMESIZE * 2 );
    AudioUnitInitialize(au);
    AudioOutputUnitStart(au);
       // Do any additional setup after loading the view, typically from a nib.
}


- (IBAction)SwitchPressed:(id)sender
{
    if(EnhancedSwitch.on)
    {
        on=1;
    }
    else
    {
        on=0;
    }
}

- (IBAction)Sir_SwitchPressed:(id)sender
{
    if(SirenSwitch.on)
    {
        sir_on=1;
    }
    else
    {
        sir_on=0;
    }
}


- (IBAction)volumesliderAction:(id)sender
{
    [[MPMusicPlayerController applicationMusicPlayer] setVolume:self.volumeslider.value];
}

@end

