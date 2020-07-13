//
//  ViewController.h
//  CoreAudio
//
//  Created by Shankar, Nikhil on 4/4/17.
//  Copyright © 2017 default. All rights reserved.
//

#import <UIKit/UIKit.h>
#import <AVFoundation/AVFoundation.h>
#import <AudioToolbox/AudioToolbox.h>
#import <AudioUnit/AudioUnit.h>
#import <MediaPlayer/MediaPlayer.h>
@interface ViewController : UIViewController{

@public
float x,y;
}
//@property (weak,nonatomic) IBOutlet UILabel *betaLabel;
@property (weak, nonatomic) IBOutlet UISwitch *EnhancedSwitch;
@property (weak, nonatomic) IBOutlet UISwitch
*SirenSwitch;
@property (weak, nonatomic) IBOutlet UISlider *volumeslider;

//@property (weak, nonatomic) IBOutlet UILabel *NormaLabel;

@property (weak, nonatomic) IBOutlet UILabel *sirenLabel;


- (IBAction)volumesliderAction:(id)sender;

//- (IBAction)buttonPressed:(id)sender;
- (IBAction)SwitchPressed:(id)sender;
- (IBAction)Sir_SwitchPressed:(id)sender;




@end

