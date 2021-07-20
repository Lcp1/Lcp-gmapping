#ifndef ACCESSTATE_H
#define ACCESSTATE_H

namespace GMapping {
    // 对应二进制000的是内部Outside,001的是外部Inside,010的是局外Allocated
enum AccessibilityState{Outside=0x0, Inside=0x1, Allocated=0x2};
};

#endif

