#pragma once
#include <QWebFrame>
#include <QPainter>

class Control : public QObject
{
    Q_OBJECT
    Q_PROPERTY(bool enabled READ getEnabled WRITE setEnabled)

private:
    QWebFrame * frame;

    bool enabled;
    bool getEnabled(){ return enabled; }
    void setEnabled(bool state){ enabled = state; }

public:
    Control(QWebFrame * frame) : frame(frame), enabled(true)
    {
        this->setObjectName("Control");
    }

    Q_INVOKABLE int callableEverywhere(){ return 1; }

public slots:
    void receiveData(QVariantMap data);

signals:
    void sendData(QVariantMap data);
};
