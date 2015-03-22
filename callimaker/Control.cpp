#include "Control.h"

void Control::receiveData(QVariantMap data)
{
    for(auto key : data.keys())
    {
        auto value = data[key];
        qDebug() << value;
    }
}
