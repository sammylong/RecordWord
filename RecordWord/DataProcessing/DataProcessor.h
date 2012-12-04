//
//  DataProcessor.h
//  AudioRecorder
//
//  Created by Sammy Long on 12/9/28.
//  Copyright (c) 2012年 Apexlearn Inc. All rights reserved.
//

#ifndef AudioRecorder_DataProcessor_h
#define AudioRecorder_DataProcessor_h

void MFCC(uint8_t *wave, size_t wave_size);
void FourierTransformTest(void);
void cosine_transform_test(void);

#endif
