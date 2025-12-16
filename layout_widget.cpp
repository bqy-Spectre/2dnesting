#include "layout_widget.h"
#include <QPainterPath>

void LayoutWidget::initialize_transform() {
    globalTransform.reset();
    auto r = rect();
    qreal scaleX = (r.width() - 5) / sheet.width();
    qreal scaleY = (r.height() - 5) / sheet.height();
    qreal scale = qMin(scaleX, scaleY);
    globalTransform.translate(r.center().rx() - sheet.center().rx() * scale,
        r.center().ry() - sheet.center().ry() * scale);
    globalTransform.scale(scale, scale);
}

void LayoutWidget::set_sheet_rect(qreal w, qreal h) {
    isCircle = false;
    sheet = QRectF(0, 0, w, h);
}

void LayoutWidget::set_sheet_circle(qreal radius) {
    isCircle = true;
    qreal d = radius * 2.0;
    sheet = QRectF(0, 0, d, d);
}

void LayoutWidget::initializeGL() {
    initializeOpenGLFunctions();
    initialize_transform();
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
}

void LayoutWidget::paintEvent(QPaintEvent* event) {
    QPainter painter;
    painter.begin(this);
    painter.setRenderHint(QPainter::Antialiasing);
    painter.fillRect(rect(), Qt::white);
    painter.save();

    qreal scaledPenWidth = 2.0 / globalTransform.m11();
    QPen scaledPen = painter.pen();
    if (scaledPenWidth == 0) {
        scaledPenWidth = 1;
    }
    scaledPen.setWidthF(scaledPenWidth);

    // 设置sheet边界线更宽
    QPen sheetPen = scaledPen;
    sheetPen.setWidthF(sheetPen.widthF() * 1.5);

    // 绘制sheet，根据类型
    painter.setTransform(globalTransform);
    painter.setPen(sheetPen);
    painter.setBrush(Qt::NoBrush);
    if (isCircle) {
        painter.drawEllipse(sheet);
    } else {
        painter.drawRect(sheet);
    }

    // Clip to boundary
    QPainterPath clipPath;
    if (isCircle) {
        clipPath.addEllipse(sheet);
    } else {
        clipPath.addRect(sheet);
    }
    painter.setClipPath(clipPath);

    // 绘制多边形
    QPen polygonPen = scaledPen;
    painter.setPen(polygonPen);
    painter.setBrush(Qt::gray);
    for (auto& p : layout) {
        painter.drawPolygon(p);
    }

    painter.restore();
    painter.end();
}

void LayoutWidget::wheelEvent(QWheelEvent* event) {
    // 处理鼠标滚轮事件来实现缩放
    // 获取滚轮滚动的角度
    int delta = event->angleDelta().y();
    QPointF mousePos = event->position();
    auto mapPt = globalTransform.inverted().map(mousePos);
    qreal mapX = mapPt.x();
    qreal mapY = mapPt.y();
    // 计算缩放因子
    qreal scaleFactor = 1.1;  // 缩放因子可以根据需要调整
    if (event->angleDelta().y() < 0) {
        scaleFactor = 1.0 / scaleFactor;  // 向下滚动进行缩小
    }
    globalTransform.translate(mapX, mapY);
    globalTransform.scale(scaleFactor, scaleFactor);
    globalTransform.translate(-mapX, -mapY);
    update();
}

void LayoutWidget::mousePressEvent(QMouseEvent* event) {
    if (event->button() == Qt::LeftButton) {
        // 记录鼠标点击位置
        m_lastPos = event->pos();
    }
}

void LayoutWidget::mouseMoveEvent(QMouseEvent* event) {
    if (event->buttons() & Qt::LeftButton) {
        // 计算鼠标移动的距离
        QPointF delta = event->pos() - m_lastPos;
        // 更新平移的距离
        globalTransform.translate(delta.x() / globalTransform.m11(),
            delta.y() / globalTransform.m22());
        // 重新绘制
        update();
        // 更新lastPos
        m_lastPos = event->pos();
    }
}

void LayoutWidget::layoutUpdate(QTableWidgetItem* n, QTableWidgetItem* o) {
    qDebug() << "layoutUpdate START";
    if (n != nullptr) {
        auto i = n->tableWidget()->item(n->row(), 0);
        auto v = i->data(Qt::UserRole).value<QList<QPolygonF>>();
        auto len = i->data(Qt::DisplayRole).value<qreal>();
        this->layout = v;
        this->length = len;
        initialize_transform();
        update();
    }
}
