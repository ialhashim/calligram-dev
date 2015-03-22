#include "WebPage.h"
#include <QDebug>

WebPage::WebPage(QObject *parent): QWebPage(parent){
    qDebug() << "Creating web page..";
}

void WebPage::javaScriptConsoleMessage(const QString& message, int lineNumber, const QString& sourceID){
   QString logEntry = message +" on line:"+ QString::number(lineNumber) +" Source:"+ sourceID;
   qDebug() << logEntry;
}
